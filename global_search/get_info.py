"""
This script is mainly for retrieve infomation enough for following analysis

"""
from global_search.thirty_party.EntrezDownloader import EntrezDownloader
from global_search.classification_script import _classificated
import random
from Bio import Entrez
from tqdm import tqdm
from Bio import Entrez, SeqIO
import io
from collections import defaultdict
import multiprocessing as mp
from bs4 import BeautifulSoup
import pandas as pd
from os.path import exists, dirname, join
import os
import click
from ete3 import NCBITaxa

ncbi = NCBITaxa()


def parse_id(infile, columns=1):
    id_list = []
    id2info = {}
    for row in tqdm(open(infile, 'r')):
        if row:
            id = row.split('\t')[columns].strip().strip('\n')
            if '||' in id:
                id = id.split('||')[-1]
            id_list.append(id)
            id2info[id] = ';'.join(
                row.split('\t')[columns+1:]).strip('\n')
    return id_list, id2info


def parse_bioproject_xml(xml_text):
    result_bucket = []
    soup = BeautifulSoup(xml_text, 'xml')
    split_out = soup.find_all("DocumentSummary")
    for each_record in split_out:
        bioproject2info = defaultdict(dict)
        body = each_record.Project
        uid = each_record['uid']
        if body is None:
            result_bucket.append(uid)
            continue
        major_parts = {_.name: _ for _ in body.children if _.name}
        PID_part = major_parts['ProjectID']
        bioproject_id = PID_part.ArchiveID['accession']
        bioproject_gid = int(PID_part.ArchiveID['id'])
        PID_des = major_parts['ProjectDescr']
        des_text = PID_des.Description.text if PID_des.Description is not None else ''
        des_title = PID_des.Title.text if PID_des.Title is not None else ''
        PID_pubmed_id = ';'.join([_['id']
                                  for _ in PID_des.find_all('Publication')])
        list_biosample = PID_des.find_all('LocusTagPrefix')
        biosample_text = ';'.join([_.get('biosample_id', 'no biosample ID')
                                   for _ in list_biosample])

        PID_type = major_parts['ProjectType']
        biological_properties = PID_type.find('BiologicalProperties')
        if biological_properties is not None:
            env_data = biological_properties.find('Environment')
            if env_data:
                for _data in env_data:
                    if _data != '\n':
                        bioproject2info[bioproject_id][_data.name] = _data.text
            Morphology_data = biological_properties.find('Morphology')
            if Morphology_data:
                for _data in Morphology_data:
                    if _data != '\n':
                        bioproject2info[bioproject_id][_data.name] = _data.text
        target_pro = PID_type.find('target')
        if target_pro:
            for _data,v in target_pro.attrs.items():
                bioproject2info[bioproject_id][_data] = v
        bioproject2info[bioproject_id]['GI'] = bioproject_gid
        bioproject2info[bioproject_id]['description'] = des_text
        bioproject2info[bioproject_id]['title'] = des_title
        bioproject2info[bioproject_id]['pubmed GI'] = PID_pubmed_id
        bioproject2info[bioproject_id]['relative biosample'] = biosample_text
        bioproject2info[bioproject_id]['number of biosamples'] = len(
            list_biosample)
        pl = each_record.ProjectLinks
        if pl is not None:
            memids_DOM = [_ for _ in pl.find_all("MemberID")]
            memberids = ';'.join([_['id'] for _ in memids_DOM])
            member_accessions = ';'.join([_['accession'] for _ in memids_DOM])
            bioproject2info[bioproject_id]['member bioproject GI'] = memberids
            bioproject2info[bioproject_id]['member bioproject accession'] = member_accessions
        result_bucket.append(bioproject2info)
    return result_bucket


def parse_biosample_xml(xml_text):
    result_bucket = []
    soup = BeautifulSoup(xml_text, 'xml')
    split_out = soup.find_all("BioSample")
    for each_record in split_out:
        biosample2info = defaultdict(dict)
        uid = each_record['id']
        accession = each_record['accession']
        biosample2info[accession]['access status'] = each_record['access']
        biosample2info[accession]['last_update'] = each_record['last_update']
        biosample2info[accession]['publication_date'] = each_record['publication_date']
        biosample2info[accession]['submission_date'] = each_record['submission_date']

        for dom in each_record.find_all('Id'):
            for attr in dom.attrs:
                if attr.startswith('db'):
                    biosample2info[accession][dom[attr]] = dom.text
        for dom in each_record.find('Description').children:
            if dom.name is not None:
                biosample2info[accession][dom.name] = dom.text
        for dom in each_record.find_all('Attribute'):
            biosample2info[accession]['attribute:' +
                                      dom['attribute_name']] = dom.text
        result_bucket.append(biosample2info)
    return result_bucket


def main(infile, odir, batch_size, fectch_size,test=False,just_seq=False):
    fectch_size = int(fectch_size)
    
    edl = EntrezDownloader(
        # An email address. You might get blocked by the NCBI without specifying one.
        email='l0404th@gmail.com',
        # An API key. You can obtain one by creating an NCBI account. Speeds things up.
        api_key='ccf9847611deebe1446b9814a356f14cde08',
        num_threads=30,                       # The number of parallel requests to make
        # The number of IDs to fetch per request
        batch_size=batch_size,
        pbar=True                             # Enables a progress bar, requires tqdm package
    )
    if not exists(odir):
        os.makedirs(odir)
    order_id_list, id2annotate = parse_id(infile, 0)
    id_list = list(set(order_id_list))
    if test:
        id_list = random.sample(id_list, 1000)

    pid2info_dict = defaultdict(dict)
    tqdm.write('from protein Accession ID to GI')
    results, failed = edl.esearch(db='protein',
                                  ids=id_list,
                                  result_func=lambda x: Entrez.read(io.StringIO(x))['IdList'])
    all_GI = list(set(results[::]))
    tqdm.write('get pid summary from each one')
    results, failed = edl.esummary(db='protein',
                                   ids=all_GI,
                                   result_func=lambda x: Entrez.read(
                                       io.StringIO(x)))
    if failed:
        tqdm.write("failed retrieve %s summary of protein ID" % len(failed))
    taxons = ['superkingdom', 'phylum', 'class',
              'order', 'family', 'genus', 'species']
    gi2pid = {}
    tqdm.write('from summary to GI and taxonomy')
    for result in tqdm(results):
        aid = result['AccessionVersion']
        pid2info_dict[aid]['GI'] = gi = result['Gi'].real
        pid2info_dict[aid]['taxid'] = result['TaxId'].real
        try:
            lineage = ncbi.get_lineage(result['TaxId'].real)
            rank = ncbi.get_rank(lineage)
            rank = {v: k for k, v in rank.items()}
            names = ncbi.get_taxid_translator(lineage)
            for c in taxons:
                if c in rank:
                    pid2info_dict[aid][c] = names[rank[c]]
        except:
            tqdm.write("failed to parse taxonomy info for ", aid)
        gi2pid[gi] = aid
    with open(join(odir,'all_gi.txt'),'w') as f1:
        f1.write('\n'.join(map(str,all_GI)))
    tqdm.write("successfully retrieve %s summary of protein ID" % len(results))
    tqdm.write('retrieving protein info')
    # if just_seq:
    #     prot_results, prot_failed = edl.efetch(db='protein',
    #                                        ids=all_GI,
    #                                        retmode='text',
    #                                        retype='fasta',
    #                                        batch_size=fectch_size,
    #                                        result_func=lambda x: list(SeqIO.parse(
    #                                            io.StringIO(x), format='fasta')))
        # with open(join(odir, 'all_seqs.faa'),'w') as f1:
        #     SeqIO.write(f1,prot_results,format='fasta-2line')
        # return 
    # else:
    prot_results, prot_failed = edl.efetch(db='protein',
                                        ids=all_GI,
                                        retmode='text',
                                        retype='gb',
                                        batch_size=fectch_size,
                                        result_func=lambda x: list(SeqIO.parse(
                                            io.StringIO(x), format='genbank')))
    if prot_failed:
        tqdm.write("failed retrieve %s genbank of protein ID" % len(failed))
    refs = ['reference_' + str(_+1) + _suffix for _ in range(10)
            for _suffix in ['',' journal',' author']]
    new_columns = ['protein accession',
                       'annotated as',
                       'seq',
                       'org',
                       'source',
                       'BioProject',
                       'BioSample',
                       'GI',
                       'taxid',
                       'nuccore id',
                       'keywords',
                       'comments'] + taxons + refs
    #all_cols = list(list(pid2info_dict.values())[0].keys())
    with open(join(odir, 'protein2INFO.tab'), 'w') as f1:
        print('\t'.join(new_columns),file=f1)
        tqdm.write('write into a dictinoary and also write into a file')
        def write_in(t):
            f1.write(t.replace('\n',' ')+'\t')
        for prot_t in tqdm(prot_results):
            aid = prot_t.id
            if aid not in id_list:
                print('error ', aid)
                if aid.split('.')[0] in id_list:
                    aid = [_ for _ in id_list if _ in aid][0]
                    pass
                else:
                    continue
            annotations = prot_t.annotations
            ref_texts = [_
                         for _ in annotations.get('references', [])
                        if 'Direct' not in _.title and _.title]
            f1.write(f'{aid}\t')
            f1.write(id2annotate.get(aid, '')+'\t')
            f1.write(str(prot_t.seq)+'\t')
            write_in(annotations.get('organism',''))
            write_in(annotations.get('source',''))
            db_ = dict([_.split(':') for _ in prot_t.dbxrefs if ':' in _])
            write_in(db_.get('BioProject',''))
            write_in(db_.get('BioSample',''))
            write_in(annotations.get('GI',''))
            write_in(annotations.get('taxid',''))
            write_in(annotations.get('db_source', '').split(' ')[-1])
            write_in(';'.join(annotations.get('keywords', [])))
            write_in(annotations.get('comment', ''))
            for t in taxons:
                f1.write(pid2info_dict.get(aid,{}).get(t,'')+'\t')
            for idx in range(10):
                if idx < len(ref_texts):
                    ref_t = ref_texts[idx]
                    write_in('\t'.join([ref_t.title,ref_t.journal,ref_t.authors]))
                else:
                    f1.write('\t'.join(['','','']))
            f1.write('\n')
            f1.flush()
            if just_seq:
                continue
            
            for idx, ref_t in list(enumerate(ref_texts))[:10]:
                pid2info_dict[aid]['reference_'+str(int(idx)+1)] = ref_t.title
                pid2info_dict[aid]['reference_' +
                                   str(int(idx)+1) + ' journal'] = ref_t.journal
                pid2info_dict[aid]['reference_' +
                                   str(int(idx)+1) + ' author'] = ref_t.authors
            pid2info_dict[aid]['nuccore id'] = annotations.get(
                'db_source', '').split(' ')[-1]
            pid2info_dict[aid]['source'] = annotations['source']
            pid2info_dict[aid]['org'] = annotations['organism']
            pid2info_dict[aid]['keywords'] = ';'.join(
                annotations.get('keywords', []))
            pid2info_dict[aid]['comments'] = annotations.get('comment', '')
            pid2info_dict[aid]['seq'] = str(prot_t.seq)
            pid2info_dict[aid].update(dict([_.split(':')
                                            for _ in prot_t.dbxrefs if ':' in _]))
            pid2info_dict[aid]['annotated as'] = [id2annotate.get(_, '')
                                                  for _ in pid2info_dict.keys()]
        if just_seq:
            tqdm.write('only perform sequence searching... completed')
            return
    if just_seq:
        tqdm.write('')
    else:
        tqdm.write(
            'processing pid to bioproject and retrieving the info of bioproject')
        set_bioprojects = list(set([_.get('BioProject', '')
                                    for _ in pid2info_dict.values()]))
        set_bioprojects = [_
                        for _ in set_bioprojects
                        if str(_) != 'nan' and _]
        results, failed = edl.esearch(db='bioproject',
                                    ids=set_bioprojects,
                                    result_func=lambda x: Entrez.read(io.StringIO(x))['IdList'])
        all_GI = results[::]
        results, failed = edl.efetch(db='bioproject',
                                    ids=all_GI,
                                    retmode='xml',
                                    retype='xml',
                                    # batch_size=1,
                                    result_func=lambda x: parse_bioproject_xml(x))
        _t = {}
        for _ in results:
            if isinstance(_, dict):
                _t.update(_)
        bioproject_df = pd.DataFrame.from_dict(_t, orient='index')
        bioproject_df.loc[:, 'GI'] = bioproject_df.loc[:, 'GI'].astype(int)
        #bioproject_df = bioproject_df.reindex(pid2info_df.loc[:, 'BioProject'])
        bioproject_df = bioproject_df.applymap(
            lambda x: x.replace('\n', ' ') if isinstance(x, str) else x)
        bioproject_df.to_excel(join(odir, 'bioproject2info.xlsx'),
                            index=1, index_label='bioproject ID')

        tqdm.write('processing pid to biosample and retrieving the info of biosample')
        set_biosamples = list(set([_.get('BioSample', '')
                                for _ in pid2info_dict.values()]))
        set_biosamples = [_
                        for _ in set_biosamples
                        if str(_) != 'nan' and _]

        results, failed = edl.esearch(db='biosample',
                                    ids=set_biosamples,
                                    result_func=lambda x: Entrez.read(io.StringIO(x))['IdList'])
        all_GI = results[::]
        results, failed = edl.efetch(db='biosample',
                                    ids=all_GI,
                                    retmode='xml',
                                    retype='xml',
                                    # batch_size=1,
                                    result_func=lambda x: parse_biosample_xml(x))
        _t = {}
        for r in results:
            _t.update(r)
        biosample_df = pd.DataFrame.from_dict(_t, orient='index')
        #biosample_df = biosample_df.reindex(pid2info_df.loc[:, 'BioSample'])
        biosample_df = biosample_df.applymap(
            lambda x: x.replace('\n', ' ') if isinstance(x, str) else x)

        biosample_df = _classificated(biosample_df)
        biosample_df.to_excel(join(odir, 'biosample2info.xlsx'),
                            index=1, index_label='biosample ID')
        # merged thems
        protein2info_df = pd.read_csv(join(odir, 'protein2INFO.tab'),sep='\t',index_col=0)
        _bioproject_df = bioproject_df.reindex(protein2info_df.loc[:,'BioProject'])
        _biosample_df = biosample_df.reindex(protein2info_df.loc[:,'BioSample'])
        _bioproject_df.index = protein2info_df.index
        _biosample_df.index = protein2info_df.index
        full_df = pd.concat([protein2info_df,_bioproject_df,_biosample_df],axis=1)
        
        full_df.to_excel(join(odir, 'full_info.xlsx'),
                            index=1, index_label='protein accession')
        
        

@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id and its annotation.')
@click.option('-o', 'odir', help='output directory')
@click.option('-bs', 'batch_size', help='number of sample fetch at each query', default=500, required=False)
@click.option('-fs', 'fectch_size', help='number of sample fetch at each query', default=50, required=False)
@click.option('-debug', 'test', help='test?', default=False, required=False, is_flag=True)
@click.option('-only_seq', 'just_seq', help='only retrieve seq?', default=False, required=False, is_flag=True)
def cli(infile, odir, test, batch_size,just_seq,fectch_size):
    batch_size = int(batch_size)
    main(infile, odir, batch_size, fectch_size,test,just_seq)


if __name__ == "__main__":
    cli()
