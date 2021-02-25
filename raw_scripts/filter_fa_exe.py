import io
import os
from os.path import *
from subprocess import check_call

import numpy as np
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
from tqdm import tqdm

from bin.ncbi_convertor import edl


def get_len(records):
    return [len([bp for bp in _ if bp != '-'])
            for _ in records]


def run(cmd):
    check_call(cmd, shell=True, stdout=open('/dev/null'))


kofam_scan = '/home-user/thliao/software/kofamscan/exec_annotation'

gene_info = {'kegg': {'nxrA': 'K00370',
                      'nxrB': 'K00371',
                      'hao': 'K10535',
                      'amoA': 'K10944',
                      'amoB': 'K10945',
                      'amoC': 'K10946'},
             'TIGFAM': {'nxrA': '',
                        'nxrB': '',
                        'hao': 'TIGR01703',
                        'amoA': 'TIGR03080',
                        'amoB': 'TIGR03079',
                        'amoC': 'TIGR03078'}}


def filter_by_relative_pos(full_df):
    pid2nuccore_ids = dict(zip(full_df.index, full_df.loc[:, 'nuccore ID']))
    nuccore_ids = [v for k, v in pid2nuccore_ids.items()]
    _nuccore_summary = []
    nuccore_summary, failed = edl.esummary(db='nuccore',
                                           ids=nuccore_ids,
                                           result_func=lambda x: Entrez.read(io.StringIO(x)))
    if failed:
        _nuccore_summary, failed = edl.esummary(db='nuccore',
                                                ids=failed,
                                                result_func=lambda x: Entrez.read(io.StringIO(x)),
                                                batch_size=1)
    nid2length = {}
    for result in nuccore_summary + _nuccore_summary:
        nid2length[result['AccessionVersion']] = result['Length'].real
    nid2length = {k: v for k, v in nid2length.items() if k in pid2nuccore_ids.values()}

    tqdm.write('%s need to manual adjust....' % len(set(pid2nuccore_ids.values()).difference(set(nid2length))))
    near_end_protein = []
    for idx, row in full_df.iterrows():
        end = row['start']
        end = row['end']
        start = row['start']
        nuccore_ID = pid2nuccore_ids[idx]
        length = nid2length.get(nuccore_ID, '')
        if length:
            if abs(end - length) <= 100 or start <= 100:
                tqdm.write('near end :' + idx)
                near_end_protein.append(idx)
        else:
            return_v = input(
                f'nuccore ID is {nuccore_ID}; does not get any info, maybe wgs, please check it manually. protein id is {idx}, end at {end},start at {start}. Y/y for indicated it near the end. ')
            if return_v.lower() == 'y':
                near_end_protein.append(idx)
            else:
                pass
    return near_end_protein


def filter_fa_by_length_dis(in_fa, ofile=None, output_records=True, down_threshold=25, upper_threshold=100, hard_filter=None):
    records = [_ for _ in SeqIO.parse(in_fa, format='fasta')]
    length_dis = [len(_.seq) for _ in records]
    down_len = np.percentile(length_dis, down_threshold)
    upper_len = np.percentile(length_dis, upper_threshold)
    print('down len: ', down_len)
    print('upper len: ', upper_len)
    _s = sorted(records, key=lambda x: len(x.seq))
    print('longest seq is %s, has %s AA' % (_s[-1].id, len(_s[-1].seq)))
    print('shortest seq is %s, has %s AA' % (_s[0].id, len(_s[0].seq)))
    if hard_filter is None:
        remained_records = [_
                            for _ in records
                            if len(_.seq) > down_len and len(_.seq) < upper_len]
    else:
        remained_records = [_
                            for _ in records
                            if len(_.seq) > hard_filter]
    print('ori number of sequences: ', len(records))
    print('remained number of sequences: ', len(remained_records))
    _s = sorted(remained_records, key=lambda x: len(x.seq))
    print('longest seq is %s, has %s AA' % (_s[-1].id, len(_s[-1].seq)))
    print('shortest seq is %s, has %s AA' % (_s[0].id, len(_s[0].seq)))
    if (ofile is not None) and (not output_records):
        if not exists(dirname(ofile)):
            os.makedirs(dirname(ofile))
        with open(ofile, 'w') as f1:

            SeqIO.write(remained_records, f1, format='fasta-2line')
    elif output_records:
        return remained_records
    return


def filter_archaea(full_df, remove_nc=False):
    # filter out archaea and NC10
    if remove_nc:
        remained_ids = full_df.index[(full_df.loc[:, 'superkingdom'] == 'Bacteria') & (~full_df.loc[:, 'phylum'].str.contains('NC10').fillna(False))]
    else:
        remained_ids = full_df.index[full_df.loc[:, 'superkingdom'] == 'Bacteria']
    return remained_ids


def remove_pralog(records, paralog_file):
    remained_records = records
    if exists(paralog_file) and paralog_file:
        all_ids = [[_ for _ in row.split(' ') if _][2] for row in open(paralog_file).read().split('\n') if not row.startswith('#') and row]
        remained_records = [_ for _ in records if _.id not in all_ids]
    return remained_records


def cluster_fa(infa, odir=None):
    if not '/' in infa:
        infa = './' + infa
    if not odir:
        odir = dirname(abspath(infa))
    if not exists(odir):
        os.makedirs(odir)
    for threshold in [90, 95, 98]:
        run(f"cd-hit -i {infa} -o {odir}/cluster_{threshold} -c 0.{threshold} -n 5  -T 20 -d 0")


def compared(infile):
    records_ori = list(SeqIO.parse(infile, format='fasta'))

    ori_medi = np.median(get_len(records_ori))
    ori_mean = np.mean(get_len(records_ori))
    ori_std = np.std(get_len(records_ori))

    print(f'ori median: {ori_medi}')
    print(f'ori mean: {ori_mean}')
    print(f'ori std: {ori_std}')


tds = ['nr_retrieve_amoB', 'nr_retrieve_amoC', 'with_genome_amoA', 'nr_retrieve_hao', 'nr_retrieve_nxrA',
       "nr_retrieve_nxrB"]
for target_dir in tds:
    g = target_dir.split('_')[-1]
    # if g == 'amoB':
    #     # target_dir = './rough_amoB'
    #     pro2full_tab = f'{target_dir}/info_dir/pro2full_info.tab'
    #     fa = f'{target_dir}/used.faa'
    if g in ['amoB', 'amoC']:
        pro2full_tab = f'{target_dir}/filtered_by_kegg.faa_aln.dir/iqtree.treefile/info_dir/pro2full_info.tab'
        fa = f'{target_dir}/filtered_by_kegg.faa'
    elif g == 'amoA':
        pro2full_tab = f'{target_dir}/used.fasta_aln.dir/iqtree.treefile/info_dir/pro2full_info.tab'
        fa = f'{target_dir}/filtered_by_kegg.faa'
    elif g == 'hao':
        pro2full_tab = f'{target_dir}/info_dir/pro2full_info.tab'
        fa = f'{target_dir}/filtered_by_kegg.faa'
    elif g == 'nxrA':
        pro2full_tab = f'{target_dir}/info_dir/pro2full_info.tab'
        fa = f'{target_dir}/filtered_by_kegg.faa'
        paralog_file = f'{target_dir}/paralog_TIGFAM.hmmscan'
    elif g == 'nirK':
        pro2full_tab = f'{target_dir}/info_dir/pro2full_info.tab'
        fa = f'{target_dir}/used.faa'
    elif g == 'nxrB':
        pro2full_tab = f'{target_dir}/info_dir/pro2full_info.tab'
        fa = f'{target_dir}/filtered_by_kegg.faa'

    full_df = pd.read_csv(pro2full_tab, sep='\t', index_col=0)

    records = list(SeqIO.parse(fa, format='fasta'))
    near_end_protein = filter_by_relative_pos(full_df)

    remained_B_ids = filter_archaea(full_df, remove_nc=False)
    remained_B_remove_nc_ids = filter_archaea(full_df, remove_nc=True)
    num_ori = len(records)

    remained_records = [_ for _ in records if _.id in remained_B_ids]
    remained_records = [_ for _ in remained_records if _.id not in near_end_protein]
    if g == 'nxrA':
        remained_records = remove_pralog(remained_records, paralog_file)


    final_fa = f'{target_dir}/with_genome_Bacteria_intact.faa'
    with open(final_fa, 'w') as f1:
        SeqIO.write(remained_records, f1, format='fasta-2line')
    # filter_fa_by_length_dis(final_fa,ofile=final_fa,hard_filter=660,output_records=False)
    if g == 'nxrA' and len(remained_records) > 1000:
        cluster_fa(final_fa,
                   f'{final_fa}_aln.dir')
        final_fa = f'{final_fa}_aln.dir/cluster_95'

    compared(final_fa)
    print('remained %s fa' % len([_ for _ in SeqIO.parse(final_fa, format='fasta')]))
    print('original %s fa' % num_ori)

    remained_records = [_ for _ in records if _.id in remained_B_remove_nc_ids]
    remained_records = [_ for _ in remained_records if _.id not in near_end_protein]
    if g == 'nxrA':
        remained_records = remove_pralog(remained_records, paralog_file)
    final_fa2 = f'{target_dir}/with_genome_Bacteria_drop_NC10_intact.faa'
    with open(final_fa2, 'w') as f1:
        SeqIO.write(remained_records, f1, format='fasta-2line')
    # filter_fa_by_length_dis(final_fa2,ofile=final_fa2,hard_filter=660,output_records=False)
    if g == 'nxrA' and len(remained_records) > 1000:
        cluster_fa(final_fa2,
                   f'{final_fa2}_aln.dir')
        final_fa2 = f'{final_fa2}_aln.dir/cluster_90'
    compared(final_fa2)
    print('remained %s fa' % len([_ for _ in SeqIO.parse(final_fa2, format='fasta')]))
    print('original %s fa' % num_ori)

    cmd = f'python3 ~/script/evolution_relative/global_search/build_tree_exe.py {final_fa}'
    check_call(cmd, shell=1)

    cmd = f'python3 ~/script/evolution_relative/global_search/build_tree_exe.py {final_fa2}'
    check_call(cmd, shell=1)
    ####### after iqtree
    cmd = f'python3 ~/script/evolution_relative/global_search/build_tree_exe.py {final_fa} .iqtree.treefile'
    check_call(cmd, shell=1)

    cmd = f'python3 ~/script/evolution_relative/global_search/build_tree_exe.py {final_fa2} .iqtree.treefile'
    check_call(cmd, shell=1)
    cmd = f"cp -r {dirname(pro2full_tab)} {final_fa}_aln.dir/iqtree.treefile"
    check_call(cmd, shell=1)
    cmd = f"cp -r {dirname(pro2full_tab)} {final_fa2}_aln.dir/iqtree.treefile"
    check_call(cmd, shell=1)
    cmd = f"python3 ~/script/evolution_relative/global_search/reannotate_tree.py {final_fa}_aln.dir/iqtree.treefile {final_fa2}_aln.dir/iqtree.treefile"
    check_call(cmd, shell=1)

if __name__ == "__main__":
    import sys

    params = sys.argv[1:]
