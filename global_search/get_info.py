"""
This script is mainly for retrieve infomation enough for following analysis

"""
from global_search.thirty_party.EntrezDownloader import EntrezDownloader
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


def main(infile, odir, batch_size, test=False):

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
    all_GI = results[::]
    tqdm.write('get pid summary from each one')
    results, failed = edl.esummary(db='protein',
                                   ids=all_GI,
                                   result_func=lambda x: Entrez.read(
                                       io.StringIO(x)))
    if failed:
        tqdm.write("failed retrieve %s summary of protein ID" % len(failed))
    taxons = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']    
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
            tqdm.write("failed to parse taxonomy info for ",aid)
        gi2pid[gi] = aid
    tqdm.write("successfully retrieve %s summary of protein ID" % len(results))
    tqdm.write('retrieving protein info')
    prot_results, prot_failed = edl.efetch(db='protein',
                                           ids=all_GI,
                                           retmode='text',
                                           retype='gb',
                                           batch_size=50,
                                           result_func=lambda x: list(SeqIO.parse(
                                               io.StringIO(x), format='genbank')))
    if prot_failed:
        tqdm.write("failed retrieve %s genbank of protein ID" % len(failed))

    with open(join(odir, 'protein2INFO.tab'),'w') as f1:
        for prot_t in prot_results:
            aid = prot_t.id
            if aid not in id_list:
                print('error ', aid)
                continue
            annotations = prot_t.annotations
            ref_texts = [_
                        for _ in annotations.get('references', [])
                        if 'Direct' not in _.title and _.title]
            for idx, ref_t in enumerate(ref_texts):
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
            pid2info_dict[aid]['annotated as'] = [id2annotate.get(_,'')
                                        for _ in pid2info_dict.keys()]

        refs = list(sorted(list(set([_ for v in pid2info_dict.values() for _ in v.keys() if _.startswith('reference')]))))
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
                    'comments'] +taxons+ refs
        f1.write('\t'.join(new_columns)+'\n')
        for aid,info_dict in pid2info_dict.items():
            print(f'{aid}\t' + '\t'.join([str(info_dict.get(_,'')).replace('\n', ' ') for _ in new_columns]),file=f1)
        
        tqdm.write('transforming dictionary into a DataFrame. ')
    pid2info_df = pd.DataFrame.from_dict(pid2info_dict, orient='index')
    pid2info_df = pid2info_df.applymap(
        lambda x: x.replace('\n', ' ') if isinstance(x, str) else x)
    pid2info_df.loc[:, 'annotated as'] = [id2annotate.get(_,'')
                                        for _ in pid2info_df.index]
    
    pid2info_df = pid2info_df.reindex(columns=new_columns)
    pid2info_df.to_excel(join(odir, 'protein2INFO.xlsx'),
                        index=1, index_label='protein accession')

    tqdm.write(
        'processing pid to bioproject and retrieving the info of bioproject')
    set_bioprojects = list(pid2info_df.loc[:,'BioProject'].unique())
    set_bioprojects = [_ for _ in set_bioprojects if str(_)!='nan']
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
    bioproject_df = bioproject_df.reindex(pid2info_df.loc[:, 'BioProject'])
    bioproject_df = bioproject_df.applymap(
        lambda x: x.replace('\n', ' ') if isinstance(x, str) else x)
    bioproject_df.to_excel(join(odir, 'bioproject2info.xlsx'),
                           index=1, index_label='bioproject ID')

    tqdm.write('processing pid to biosample and retrieving the info of biosample')
    set_biosamples = [_ for _ in list(
        pid2info_df.loc[:, 'BioSample'].unique()) if isinstance(_, str)]

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
    biosample_df = biosample_df.reindex(pid2info_df.loc[:, 'BioSample'])
    biosample_df = biosample_df.applymap(
        lambda x: x.replace('\n', ' ') if isinstance(x, str) else x)
    biosample_df.to_excel(join(odir, 'biosample2info.xlsx'),
                          index=1, index_label='biosample ID')


@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id and its annotation.')
@click.option('-o', 'odir', help='output directory')
@click.option('-bs', 'batch_size', help='number of sample fetch at each query', default=500, required=False)
@click.option('-debug', 'test', help='test?', default=False, required=False, is_flag=True)
def cli(infile, odir, test, batch_size):
    batch_size = int(batch_size)
    main(infile, odir, batch_size,test)


if __name__ == "__main__":
    cli()
