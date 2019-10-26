"""
This script is mainly for retrieve infomation enough for following analysis

"""
from global_search.thirty_party.EntrezDownloader import EntrezDownloader
from global_search.classification_script import _classificated
from global_search.thirty_party.metadata_parser import *
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
import pickle

ncbi = NCBITaxa()

taxons = ['superkingdom', 'phylum', 'class',
          'order', 'family', 'genus', 'species']
source_file = '/home-user/thliao/resource/NCBI2habitat.csv'
tmp_dir = './.tmp_getINFO'
def access_intermedia(obj=None,ofile=None):
    if (ofile is not None) and (obj is not None):
        if not exists(dirname(ofile)):
            os.makedirs(dirname(ofile))
        pickle.dump(obj,open(ofile,'wb'))
    elif (not obj) and ofile:
        tqdm.write("Detecting cache file. read it")
        return pickle.load(open(ofile,'rb'))
    
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


def get_Normal_ID(id_list, fetch_size=30, edl=None):
    _md5 = str(hash(';'.join(sorted(id_list))))
    if exists(join(tmp_dir,_md5)):
        id2gi = access_intermedia(ofile=join(tmp_dir,_md5))
    else:
        tqdm.write('from protein Accession ID to GI')
        results, failed = edl.esearch(db='protein',
                                    ids=id_list,
                                    result_func=lambda x: Entrez.read(io.StringIO(x))['IdList'],
                                    batch_size=1
                                    )
        id2gi = dict(results)
        access_intermedia(id2gi,ofile=join(tmp_dir,_md5))
    gi2id = {v:k for k,v in id2gi.items()}
    all_GI = list(set(id2gi.values()))
    tqdm.write('get pid summary from each one')
    results, failed = edl.esummary(db='protein',
                                   ids=all_GI,
                                   result_func=lambda x: Entrez.read(
                                       io.StringIO(x)))
    if failed:
        tqdm.write("failed retrieve %s summary of protein ID" % len(failed))
    tqdm.write('from summary to GI and taxonomy')
    pid2info_dict = defaultdict(dict)
    for result in tqdm(results):
        aid = result['AccessionVersion']
        if aid not in gi2id:
            _aid = [_ for _ in id2gi if _ in aid]
            if not _aid:
                continue
                #print('error from esummary with gi', aid)
            else:
                right_aid = _aid[0]
        else:
            right_aid = aid
        pid2info_dict[right_aid]['GI'] = result['Gi'].real
        pid2info_dict[right_aid]['taxid'] = result['TaxId'].real
        pid2info_dict[right_aid]['accession'] = aid
        try:
            lineage = ncbi.get_lineage(result['TaxId'].real)
            rank = ncbi.get_rank(lineage)
            rank = {v: k for k, v in rank.items()}
            names = ncbi.get_taxid_translator(lineage)
            for c in taxons:
                if c in rank:
                    pid2info_dict[right_aid][c] = names[rank[c]]
        except:
            tqdm.write("failed to parse taxonomy info for ", aid)
            
    tqdm.write("successfully retrieve %s summary of protein ID" % len(results))
    tqdm.write('retrieving protein info')
    if exists(join(tmp_dir,_md5+'_efetch')):
        pid2gb = access_intermedia(ofile=join(tmp_dir,_md5+'_efetch'))
    else:
        prot_results, prot_failed = edl.efetch(db='protein',
                                            ids=all_GI,
                                            retmode='text',
                                            retype='gb',
                                            batch_size=fetch_size,
                                            result_func=lambda x: list(SeqIO.parse(
                                                io.StringIO(x), format='genbank')))
        if prot_failed:
            tqdm.write("failed retrieve %s genbank of protein ID" % len(prot_failed))
        pid2gb = {}
        for record in prot_results:
            right_aid = [k for k,v in pid2info_dict.items() if v['accession']==record.id]
            if not right_aid:
                print('Wrong ????',record.id)
            else:
                pid2gb[right_aid[0]] = record
        if set(pid2info_dict).difference(pid2gb):
            remained_id = list(set(pid2info_dict).difference(pid2gb))
            remained_id_gis = [id2gi[_] for _ in remained_id]
            prot_results, prot_failed = edl.efetch(db='protein',
                                            ids=remained_id_gis,
                                            retmode='text',
                                            retype='gb',
                                            batch_size=1,
                                            result_func=lambda x: list(SeqIO.read(
                                                io.StringIO(x), format='genbank')))
            for record in prot_results:
                right_aid = [k for k,v in pid2info_dict.items() if v['accession']==record.id]
                if not right_aid:
                    print('Wrong ????',record.id)
                else:
                    pid2gb[right_aid[0]] = record
        access_intermedia(pid2gb,ofile=join(tmp_dir,_md5+'_efetch'))
    return pid2gb, pid2info_dict

def get_WP_info(id_list, edl):
    def _parse_wp(t):
        whole_df = pd.read_csv(io.StringIO(t), sep='\t', header=None)
        aid = whole_df.iloc[0, 6]
        return [(aid, whole_df)]
    _md5 = str(hash(';'.join(sorted(id_list))))
    if exists(join(tmp_dir,_md5)):
        id2gi = access_intermedia(ofile=join(tmp_dir,_md5))
    else:
        tqdm.write('from protein Accession ID to GI')
        results, failed = edl.esearch(db='protein',
                                    ids=id_list,
                                    result_func=lambda x: Entrez.read(io.StringIO(x))['IdList'],
                                    batch_size=1)
        id2gi = dict(results)
        access_intermedia(id2gi,ofile=join(tmp_dir,_md5))
    gi2id = {v:k for k,v in id2gi.items()}
    all_GI = list(set(id2gi.values()))

    tqdm.write('get pid summary from each one')
    results, failed = edl.efetch(db='protein',
                                    ids=all_GI,
                                    retmode='ipg',
                                    retype='xml',
                                    batch_size=1,
                                    result_func=lambda x: _parse_wp(x))
    failed_id = []
    aid2info = {}
    for (aid, aid_df) in results:
        if aid not in id2gi:
            _c = [_ for _ in id2gi if aid.split('.')[0] in _]
            if not _c:
                print('unexpected ',aid)
                continue
            else:
                aid = _c[0]
        assembly_id = [_ for _ in aid_df.iloc[:, -1] if not pd.isna(_)]
        if assembly_id:
            aid2info[aid] = assembly_id[-1]
        else:
            aid2info[aid] = ''
            failed_id.append(aid)
    if failed_id:
        print("failed id %s don't have assembly id" % ';'.join(failed_id))
    assembly_id_list = list([_ for _ in aid2info.values() if _])
    results, failed = edl.esearch(db='assembly',
                                  ids=assembly_id_list,
                                  result_func=lambda x: Entrez.read(io.StringIO(x))['IdList'])
    all_GI = results[::]
    results, failed = edl.esummary(db='assembly',
                                   ids=all_GI,
                                   # batch_size=1,
                                   result_func=lambda x: parse_assembly_xml(x))
    _t = {}
    for _ in results:
        aid = list(_.keys())[0]
        v = list(_.values())[0]
        if aid not in aid2info:
            aid_s = set([_
                         for _ in aid2info.values()
                         if aid.split('_')[1].split('.')[0] in _])
            for aid in aid_s:
                _t.update({aid: v})
        else:
            _t.update({aid: v})
    pid2info = {pid: _t.get(assid, {})
                for pid, assid in aid2info.items()
                }
    return pid2info


def main(infile, odir, batch_size, fetch_size, test=False, just_seq=False, edl=None):
    fetch_size = int(fetch_size)

    if not exists(odir):
        os.makedirs(odir)
    order_id_list, id2annotate = parse_id(infile, 0)
    id_list = list(set(order_id_list))
    if test:
        id_list = random.sample(id_list, 1000)
    # for WP (refseq)
    tqdm.write('first retrieve WP/refsep protein accession. ')
    WPid_list = [_ for _ in id_list if _.startswith('WP_')]
    WPid2info = get_WP_info(WPid_list, edl=edl)

    # for others
    tqdm.write('then retrieve other protein accession. ')
    id_list = [_ for _ in id_list if not _.startswith('WP_')]

    pid2gb, pid2info_dict = get_Normal_ID(id_list, fetch_size=fetch_size, edl=edl)
    not_get_id = set(id_list).difference(set(pid2gb))
    tqdm.write(len(not_get_id),' failed retrieving...')
    tqdm.write('including ',';'.join(not_get_id))
    # init header
    refs = ['reference_' + str(_+1) + _suffix
            for _ in range(10)
            for _suffix in ['', ' journal', ' author']]
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
    pid2bioproject = {}
    pid2biosample = {}
    with open(join(odir, 'protein2INFO.tab'), 'w') as f1:
        print('\t'.join(new_columns), file=f1)
        tqdm.write('write into a dictinoary and also write into a file')

        def write_in(t):
            f1.write(t.replace('\n', ' ').replace('\t', ' ')+'\t')
        for aid in tqdm(id_list):
            if aid in pid2gb:
                prot_t = pid2gb[aid]
            else:
                _c = [_ for _ in pid2gb if _.split('.')[0] in aid]
                if not _c:
                    print('not found',aid)
                    continue
                else:
                    aid = _c[0]
                    prot_t = pid2gb[aid]
        # for prot_t in tqdm(prot_results):
            # aid = prot_t.id
            # if aid not in id_list:
            #     if aid.split('.')[0] in id_list:
            #         aid = [_ for _ in id_list if _ in aid][0]
            #         pass
            #     else:
            #         print('error ', aid)
            #         f1.write('\t'.join([aid]+['']*len(new_columns)))
            #         continue
            annotations = prot_t.annotations
            ref_texts = [_
                         for _ in annotations.get('references', [])
                         if 'Direct' not in _.title and _.title]
            f1.write(f'{aid}\t')
            f1.write(id2annotate.get(aid, '')+'\t')
            f1.write(str(prot_t.seq)+'\t')
            write_in(annotations.get('organism', ''))
            write_in(annotations.get('source', ''))
            db_ = dict([_.split(':') for _ in prot_t.dbxrefs if ':' in _])
            write_in(db_.get('BioProject', ''))
            write_in(db_.get('BioSample', ''))
            write_in(str(pid2info_dict.get(aid, {}).get('GI', '')))
            write_in(str(pid2info_dict.get(aid, {}).get('taxid', '')))
            write_in(annotations.get('db_source', '').split(' ')[-1])
            write_in(';'.join(annotations.get('keywords', [])))
            write_in(annotations.get('comment', ''))
            for t in taxons:
                f1.write(pid2info_dict.get(aid, {}).get(t, '')+'\t')
            for idx in range(10):
                if idx < len(ref_texts):
                    ref_t = ref_texts[idx]
                    write_in(
                        '\t'.join([ref_t.title, ref_t.journal, ref_t.authors]))
                else:
                    f1.write('\t'.join(['', '', '']))
            f1.write('\n')
            f1.flush()
            pid2bioproject[aid] = db_.get('BioProject', '')
            pid2biosample[aid] = db_.get('BioSample', '')

        for pid, info_dict in WPid2info.items():
            f1.write('\t'.join([pid,
                                '',
                                '',
                                info_dict.get('SpeciesName', ''),
                                '',
                                info_dict.get('BioprojectAccn', ''),
                                info_dict.get('BioSampleAccn', ''),
                                '',
                                info_dict.get('SpeciesTaxid', ''),
                                '',
                                '',
                                ''
                                ] + [''] * 30))
            pid2bioproject[pid] = info_dict.get('BioprojectAccn', '')
            pid2biosample[pid] = info_dict.get('BioSampleAccn', '')
            f1.write('\n')
            f1.flush()
    if just_seq:
        tqdm.write('only perform sequence searching... completed')
        return
    else:
        tqdm.write(
            'processing pid to bioproject and retrieving the info of bioproject')
        set_bioprojects = list(set(pid2bioproject.values()))
        set_bioprojects = [_
                           for _ in set_bioprojects
                           if _]
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
        bioproject_df = bioproject_df.applymap(
            lambda x: x.replace('\n', ' ') if isinstance(x, str) else x)
        bioproject_df.to_excel(join(odir, 'bioproject2info.xlsx'),
                               index=1, index_label='bioproject ID')

        tqdm.write(
            'processing pid to biosample and retrieving the info of biosample')
        ###
        set_biosamples = list(set(pid2biosample.values()))
        set_biosamples = [_
                          for _ in set_biosamples
                          if _]

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
        biosample_df = biosample_df.applymap(
            lambda x: x.replace('\n', ' ') if isinstance(x, str) else x)

        biosample_df = _classificated(biosample_df)
        biosample_df.to_excel(join(odir, 'biosample2info.xlsx'),
                              index=1, index_label='biosample ID')
        # merged thems
        protein2info_df = pd.read_csv(
            join(odir, 'protein2INFO.tab'), sep='\t', index_col=0)
        _bioproject_df = bioproject_df.reindex(
            protein2info_df.loc[:, 'BioProject'])
        _biosample_df = biosample_df.reindex(
            protein2info_df.loc[:, 'BioSample'])
        _bioproject_df.index = protein2info_df.index
        _biosample_df.index = protein2info_df.index
        full_df = pd.concat(
            [protein2info_df, _bioproject_df, _biosample_df], axis=1)

        if exists(source_file):
            source_df = pd.read_csv(source_file, encoding='GBK')
            source_df.loc[:, 'tmp'] = source_df.iloc[:, 0] + \
                ';'+source_df.iloc[:, 1]
            source_df = source_df.set_index('tmp')
            source_df = source_df.loc[~source_df.index.duplicated(), :]
            # bioproject
            tmp = [';'.join(list(map(str, row))[:2])
                   for row in full_df.loc[:, ['BioProject', 'BioSample']].values]

            _d1 = source_df.reindex(tmp)
            for idx, (_, v) in enumerate(full_df.iterrows()):
                if not pd.isna(_d1.iloc[idx, 2]) and pd.isna(v['habitat']):
                    # print(_,_d1.iloc[idx,2])
                    full_df.loc[_, 'habitat'] = _d1.iloc[idx, 2]
        full_df.to_excel(join(odir, 'full_info.xlsx'),
                         index=1, index_label='protein accession')


@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id and its annotation.')
@click.option('-o', 'odir', help='output directory')
@click.option('-bs', 'batch_size', help='number of sample fetch at each query', default=500, required=False)
@click.option('-fs', 'fetch_size', help='number of sample fetch at each query', default=50, required=False)
@click.option('-debug', 'test', help='test?', default=False, required=False, is_flag=True)
@click.option('-only_seq', 'just_seq', help='only retrieve seq?', default=False, required=False, is_flag=True)
def cli(infile, odir, test, batch_size, just_seq, fetch_size):
    batch_size = int(batch_size)

    edl = EntrezDownloader(
        # An email address. You might get blocked by the NCBI without specifying one.
        email='l0404th@gmail.com',
        # An API key. You can obtain one by creating an NCBI account. Speeds things up.
        api_key='ccf9847611deebe1446b9814a356f14cde08',
        num_threads=30,   # The number of parallel requests to make
        # The number of IDs to fetch per request
        batch_size=batch_size,
        pbar=True  # Enables a progress bar, requires tqdm package
    )
    main(infile, odir, batch_size, fetch_size, test, just_seq, edl=edl)


if __name__ == "__main__":
    cli()
