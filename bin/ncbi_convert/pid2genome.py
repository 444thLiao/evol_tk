"""
This script is mainly for implementing genome ID with its protein accession id.
It will generate a tab separated table.
"""
from bin.ncbi_convert import edl, access_intermedia, parse_id,unpack_gb
from bin.ncbi_convert.pid2GI import pid2GI
from bin.ncbi_convert.pid2tax import GI2tax

from os.path import exists, join, dirname
from tqdm import tqdm
from Bio import Entrez,SeqIO
import io
import os
import click
from pid2GI import pid2GI
import pandas as pd


    
    
def get_normal_ID_gb(pid2info_dict,fetch_size):
    """
    No GI info inside the genbank file retrieved... 
    with pid2info_dict which contained both exact accession ID and manual ID
    
    
    """
    all_GI = [_.get('GI','') for k,_ in pid2info_dict.items()]
    all_GI = [_ for _ in all_GI if _]
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
            tqdm.write('get unexpected record: ' + record.id + '\n')
        else:
            pid2gb[right_aid[0]] = record
    remained_pid = set(pid2info_dict).difference(set(pid2gb))
    if remained_pid:
        tqdm.write(str(len(remained_pid)) + ' proteins are missing, retrieve it one by one again. ')
        remained_GI = [pid2info_dict[_]['GI'] for _ in remained_pid]
        prot_results, prot_failed = edl.efetch(db='protein',
                                        ids=remained_GI,
                                        retmode='text',
                                        retype='gb',
                                        batch_size=1,
                                        result_func=lambda x: list(SeqIO.read(
                                            io.StringIO(x), format='genbank')))
        for record in prot_results:
            right_aid = [k for k,v in pid2info_dict.items() if v['accession']==record.id]
            if not right_aid:
                tqdm.write('get unexpected record: '+record.id)
            else:
                pid2gb[right_aid[0]] = record
    
    final_pid2gb = {}
    for pid,info_dict in pid2info_dict.items():
        info_dict = info_dict.copy()
        info_dict.update(unpack_gb(pid2gb.get(pid,{})))
        final_pid2gb[pid] = info_dict
    return final_pid2gb


def get_WP_assembly(pid2info_dict, ):
    def _parse_wp(t):
        whole_df = pd.read_csv(io.StringIO(t), sep='\t', header=None)
        aid = whole_df.iloc[0, 6]
        return [(aid, whole_df)]
    
    all_GI = [_.get('GI','') for k,_ in pid2info_dict.items()]
    all_GI = [_ for _ in all_GI if _]

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
        if aid not in pid2info_dict:
            _c = [pid for pid,v in pid2info_dict if v['accession'] == aid]
            if not _c:
                tqdm.write('get unexpected record: ' + aid)
                tqdm.write('pass it')
                continue
            else:
                aid = _c[0]
        assembly_id = [_ for _ in aid_df.iloc[:, -1] if not pd.isna(_)]
        # from last columns, get not nan one.
        if assembly_id:
            aid2info[aid] = assembly_id[-1]
        else:
            aid2info[aid] = ''
            failed_id.append(aid)
    if failed_id:
        tqdm.write("failed id %s don't have any assembly id" % ';'.join(failed_id))
    
    pid2assembly = {}
    for pid,info_dict in pid2info_dict.items():
        info_dict = info_dict.copy()
        info_dict.update({'assembly ID':aid2info.get(pid,'')} )
    return pid2assembly

def pid2genome_assembly(pid2gi,fetch_size):
    """
    could not use elink to directly from protein accession id -> assembly/biosample (So it is complicated and slowly)
    1. separate two kinds of protein ID.
        a. startwith WP, which mean refseq
        b. others
    
    """
    suffix = 'pid2genome_info'
    pid_list = list(pid2gi)
    pid2info_dict = GI2tax(pid2gi)
    _cache = access_intermedia(pid_list, suffix=suffix)
    # this is important, because it contains exact accession ID and manual provided ID.
    normal_IDs = [pid for pid in pid_list if pid and not pid.startwith('WP_')]
    refseq_IDs = [pid for pid in pid_list if pid and pid.startwith('WP_')]
    
    # normal
    normal_pid2info_dict = {pid:v for pid,v in pid2info_dict.items() if pid in normal_IDs}
    pid2gb = get_normal_ID_gb(normal_pid2info_dict,fetch_size=fetch_size)
    
    # refseq
    refseq_pid2info_dict = {pid:v for pid,v in pid2info_dict.items() if pid in refseq_IDs}
    pid2assembly = get_WP_assembly(refseq_pid2info_dict)
    
    # summarized them
    


def main(infile, ofile, force=False):
    order_id_list, id2annotate = parse_id(infile)
    id2gi = {}
    if isinstance(id2annotate[order_id_list[0]], dict):
        # it is a dict, so it contains other infomation or implemented GI. it may be passed over.
        if 'GI' in id2annotate[order_id_list[0]]:
            print("provided file already contains `GI` column(doesn't check the validation/completeness). Giving `force` param to overwrite/implement it. ")
            if not force:
                id2gi = {k:id2annotate[k]['GI'] for k in order_id_list}
        # todo: re-implemented original infomation into `ofile` from `infile`
    else:
        # no header, just a list of IDs
        pass
    if not id2gi:
        id2gi = pid2GI(order_id_list)
    
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))

    with open(ofile, 'w') as f1:
        print('#accession ID\tGI', file=f1)
        for id, GI in id2gi.items():
            print(f'{id}\t{GI}', file=f1)



@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id ')
@click.option('-o', 'ofile', help='output file')
@click.option('-f', 'force', help='force overwrite?', default=False, required=False, is_flag=True)
def cli(infile, ofile, force):
    main(infile, ofile, force)


if __name__ == "__main__":
    cli()
