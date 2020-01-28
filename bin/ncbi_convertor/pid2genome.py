"""
This script is mainly for implementing genome ID and its assembly ID with its protein accession id.
It will generate a tab separated table.
Because each protein ID would corresponding a lots of identical protein groups (which this script mainly used database)
it will generate more line than input ID. but the set of/deduplicated column 1 must match the input IDs. If not ,it will leave it blank row.

For filtration amplicons, it could simply find proteins which only contain NONE assembly ID columns. (it simple!)
header is ['accession ID', 'GI', 'assembly_ID', 'nuccore ID', 'start', 'end', 'strand']
or look like accession ID\tGI\tassembly_ID\tnuccore ID\tstart\tend\tstrand
"""
import io
import os
from os.path import exists, dirname

import click
import pandas as pd
from tqdm import tqdm

from bin.ncbi_convert import edl, access_intermedia, parse_id
from bin.ncbi_convert.pid2GI import pid2GI
from bin.ncbi_convert.pid2tax import GI2tax


def get_protein_pos_assembly_INFO(pid2info_dict, suffix='pid2genome_info'):
    def _parse_ipg(t):
        bucket = []
        whole_df = pd.read_csv(io.StringIO(t), sep='\t', header=None)
        gb = whole_df.groupby(0)
        all_dfs = [gb.get_group(x) for x in gb.groups]
        for indivi_df in all_dfs:
            indivi_df.index = range(indivi_df.shape[0])
            # aid = indivi_df.iloc[0, 6]
            indivi_df = indivi_df.fillna('')
            pos = [(row[2], row[3], row[4], row[5])
                   for row in indivi_df.values]
            gb = [row[10]
                  for row in indivi_df.values]
            for row in indivi_df.values:
                bucket.append((row[6], pos, gb))
        return bucket

    all_GI = [_.get('GI', '') for k, _ in pid2info_dict.items()]
    all_GI = [_ for _ in all_GI]
    tqdm.write('get pid summary from Identical Protein Groups summary')
    results, failed = edl.efetch(db='protein',
                                 ids=all_GI,
                                 retmode='ipg',
                                 retype='xml',
                                 result_func=lambda x: _parse_ipg(x))
    pid2assembly_dict = {}
    for pid, nuc_info, assembly_info in tqdm(results):
        pid2assembly_dict[pid] = dict(zip(assembly_info, nuc_info))
        # for nuccore which no assembly ID, it will drop it by accident. by for now, it is ok.

    pid2assembly_dict = {pid: pid2assembly_dict.get(pid2info_dict[pid].get('accession', ''), {})
                         for pid in pid2info_dict}
    assert len(pid2assembly_dict) == len(pid2info_dict)
    access_intermedia(pid2assembly_dict, suffix=suffix)
    return pid2assembly_dict


def pid2genome_assembly(pid2gi, redo=False):
    """
    could not use elink to directly from protein accession id -> assembly/biosample (So it is complicated and slowly)
    """
    suffix = 'pid2genome_info'
    pid_list = list(pid2gi)
    # get itself cache? if have, just finish
    _cache = access_intermedia(pid_list, suffix=suffix, redo=redo)
    if _cache is not None:
        return _cache
    # get pre-requested cache? if have, just go on
    _cache = access_intermedia(pid_list, suffix='pid2tax', redo=redo)
    if _cache is not None:
        pid2info_dict = _cache
    else:
        pid2info_dict = GI2tax(pid2gi, redo=redo)

    pid2assembly_nuc_info = get_protein_pos_assembly_INFO(pid2info_dict, suffix=suffix)
    return pid2assembly_nuc_info


def main(infile, ofile, force=False, redo=False):
    order_id_list, id2annotate = parse_id(infile)
    id2gi = {}
    if isinstance(id2annotate[order_id_list[0]], dict):
        # it is a dict, so it contains other infomation or implemented GI. it may be passed over.
        if 'GI' in id2annotate[order_id_list[0]]:
            print("provided file already contains `GI` column(doesn't check the validation/completeness). Giving `force` param to overwrite/implement it. ")
            if not force:
                id2gi = {k: id2annotate[k]['GI'] for k in order_id_list}
        # todo: re-implemented original infomation into `ofile` from `infile`
    else:
        # no header, just a list of IDs
        pass
    if not id2gi:
        id2gi = pid2GI(order_id_list, redo=redo)
    pid2assembly_dict = pid2genome_assembly(id2gi, redo=redo)

    if not exists(dirname(ofile)) and dirname(ofile):
        os.makedirs(dirname(ofile))

    if exists(ofile) and not force:
        tqdm.write("detect existing " + ofile + ' no force param input, so it quit instead of writing.')
        return

    with open(ofile, 'w') as f1:
        print('#accession ID\tGI\tassembly_ID\tnuccore ID\tstart\tend\tstrand', file=f1)
        for pid, nuc_dict in pid2assembly_dict.items():
            GI = id2gi[pid]
            for assembly_id, info in nuc_dict.items():
                print(f'{pid}\t{GI}\t{assembly_id}\t' + '\t'.join(map(str, info)), file=f1)
    tqdm.write('finish writing into ' + ofile)


@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id ')
@click.option('-o', 'ofile', help='output file')
@click.option('-f', 'force', help='force overwrite?', default=False, required=False, is_flag=True)
@click.option('-redo', 'redo', help='use cache or not? default is use the cache.', default=False, required=False, is_flag=True)
def cli(infile, ofile, force, redo):
    main(infile, ofile, force, redo)


if __name__ == "__main__":
    cli()

# def get_normal_ID_gb(pid2info_dict,fetch_size):
#     """
#     No GI info inside the genbank file retrieved... 
#     with pid2info_dict which contained both exact accession ID and manual ID


#     """
#     all_GI = [_.get('GI','') for k,_ in pid2info_dict.items()]
#     all_GI = [_ for _ in all_GI if _]
#     prot_results, prot_failed = edl.efetch(db='protein',
#                                         ids=all_GI,
#                                         retmode='text',
#                                         retype='gb',
#                                         batch_size=fetch_size,
#                                         result_func=lambda x: list(SeqIO.parse(
#                                             io.StringIO(x), format='genbank')))
#     if prot_failed:
#         tqdm.write("failed retrieve %s genbank of protein ID" % len(prot_failed))
#     pid2gb = {}
#     for record in prot_results:
#         right_aid = [k 
#                      for k,v in pid2info_dict.items() 
#                      if v['accession']==record.id]
#         if not right_aid:
#             tqdm.write('get unexpected record: ' + record.id + '\n')
#         else:
#             pid2gb[right_aid[0]] = record
#     remained_pid = set(pid2info_dict).difference(set(pid2gb))
#     if remained_pid:
#         tqdm.write(str(len(remained_pid)) + ' proteins are missing, retrieve it one by one again. ')
#         remained_GI = [pid2info_dict[_]['GI'] for _ in remained_pid]
#         prot_results, prot_failed = edl.efetch(db='protein',
#                                         ids=remained_GI,
#                                         retmode='text',
#                                         retype='gb',
#                                         batch_size=1,
#                                         result_func=lambda x: list(SeqIO.read(
#                                             io.StringIO(x), format='genbank')))
#         for record in prot_results:
#             right_aid = [k for k,v in pid2info_dict.items() if v['accession']==record.id]
#             if not right_aid:
#                 tqdm.write('get unexpected record: '+record.id)
#             else:
#                 pid2gb[right_aid[0]] = record

#     final_pid2gb = {}
#     tqdm.write('writing genbank into a dict, it may take much memory of your computer. If you are retrieving over 10K pid, be careful.')
#     for pid,info_dict in pid2info_dict.items():
#         info_dict = info_dict.copy()
#         info_dict.update(unpack_gb(pid2gb.get(pid,{})))
#         final_pid2gb[pid] = info_dict
#     return final_pid2gb


# def get_WP_assembly(pid2info_dict, ):
#     def _parse_wp(t):
#         whole_df = pd.read_csv(io.StringIO(t), sep='\t', header=None)
#         aid = whole_df.iloc[0, 6]
#         return [(aid, whole_df)]

#     all_GI = [_.get('GI','') for k,_ in pid2info_dict.items()]
#     all_GI = [_ for _ in all_GI if _]

#     tqdm.write('get WP/refseq pid summary one by one')
#     results, failed = edl.efetch(db='protein',
#                                     ids=all_GI,
#                                     retmode='ipg',
#                                     retype='xml',
#                                     batch_size=1,
#                                     result_func=lambda x: _parse_wp(x))
#     if failed:
#         tqdm.write('%s GI failed' % len(failed))
#     failed_id = []
#     aid2info = {}
#     for (aid, aid_df) in results:
#         if aid not in pid2info_dict:
#             _c = [pid for pid,v in pid2info_dict if v['accession'] == aid]
#             if not _c:
#                 tqdm.write('get unexpected record: ' + aid)
#                 tqdm.write('pass it')
#                 continue
#             else:
#                 aid = _c[0]
#         pinfo = {}
#         for _,row in aid_df.iterrows():
#             pinfo['assembly'] = row[10]
#             pinfo['nuccore id'] = row[2]
#             pinfo['nuc start'] = row[3]
#             pinfo['nuc end'] = row[4]
#             pinfo['nuc strand'] = row[5]
#         assembly_id = [_ for _ in aid_df.iloc[:, -1] if not pd.isna(_)]
#         # from last columns, get not nan one.
#         if assembly_id:
#             aid2info[aid] = assembly_id[-1]
#         else:
#             aid2info[aid] = ''
#             failed_id.append(aid)
#     if failed_id:
#         tqdm.write("failed id %s don't have any assembly id" % ';'.join(failed_id))

#     pid2assembly = {}
#     for pid,info_dict in pid2info_dict.items():
#         info_dict = info_dict.copy()
#         info_dict.update({'assembly ID':aid2info.get(pid,'')} )
#         pid2assembly[pid] = info_dict
#     return pid2assembly
