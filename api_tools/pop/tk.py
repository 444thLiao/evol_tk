
import copy
import io
import itertools
import os
from collections import Counter, defaultdict
from glob import glob
from os.path import *
from subprocess import check_call

import numpy as np
import pandas as pd
import plotly
import plotly.graph_objects as go
from api_tools import *
from Bio import AlignIO, SeqIO
from ete3 import Tree
from IPython.display import Image
from numpy.core.shape_base import block
from plotly.subplots import make_subplots
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm


# def fix_block_df(block_df,contig2size):
#     # since the strand would significant affect the coordinates transformation
#     # we need to fix it in advance
#     def _tmp(row, idx, end_col, start_col, strand_col, genome_col):
#         if row[strand_col] == -1:
#             c = row[genome_col].split('.')[-1]
#             if c in contig2size:
#                 size = contig2size[c]
#                 forward_idx = np.arange(size+1)
#                 reverse_idx = forward_idx[::-1]
#                 block_df.loc[idx,
#                              f'corrected_{start_col}'] = reverse_idx[row[start_col]]
#                 block_df.loc[idx,
#                              f'corrected_{end_col}'] = reverse_idx[row[end_col]]
#         else:
#             block_df.loc[idx, f'corrected_{start_col}'] = row[start_col]
#             block_df.loc[idx, f'corrected_{end_col}'] = row[end_col]
#     for idx, row in block_df.iterrows():
#         end_col, start_col, strand_col, genome_col = [
#             'End', 'Start', 'strand', 'Genome1']
#         _tmp(row, idx, end_col, start_col, strand_col, genome_col)
#         end_col, start_col, strand_col, genome_col = [
#             'End.1', 'Start.1', 'strand.1', 'Genome2']
#         _tmp(row, idx, end_col, start_col, strand_col, genome_col)
#     return block_df

def parse_sline(row,contig2length):
    # for mugsy
    # the position of the reverse complemented sequence should be corrected.
    rows = row.split(' ')
    if '\t\t' not in row:
        name = rows[1]
        start = int(rows[2])
        size = int(rows[3])
    else:
        name, start = rows[1].split('\t\t')
        start = int(start)
        size = int(rows[2])
    contig = name.rsplit('.')[-1]
    length = contig2length[contig]
    strand = 1 if rows[3] == '+' else -1
    if strand == - 1:
        start = length-start
    end = start + size
    seq = rows[-1]
    return [name, start, end, strand, seq]


def parse_maf(maf,contig2length):
    block_num = 0
    block_dict = {}
    for row in open(maf):
        row = row.strip()
        if row.startswith('a '):
            block_num += 1
            block_name = f"block{block_num}"
            block_dict[block_name] = []
        if row.startswith('s '):
            # this is a block
            block_dict[block_name].append(parse_sline(row,contig2length))
    return block_dict

def write_block(block_dict, ofile, odir, min_length=10000,min_genomes=100):
    # only block larger than min_length would be output
    if not exists(odir):
        os.system(f"mkdir -p {odir}")
    f_tab = open(ofile, 'w')
    text = f"Block\tGenome1\tStart\tEnd\tstrand\tGenome2\tStart\tEnd\tstrand\tNumber of diff bases\tSize\n"
    f_tab.write(text)
    for block, info in block_dict.items():
        if len(info) < min_genomes:
            continue
        _1 = info[0]
        if _1[2]-_1[1] < min_length:
            continue
        # 0 is name, normally is Genome.Contig_name
        # 1 is start. its original genomic coordinates
        # 2 is end. its original genomic coordinates
        for i1, i2 in itertools.combinations(info, 2):
            seq1, seq2 = i1[-1], i2[-1]
            g1, g2 = i1[0], i2[0]
            s1, e1, s2, e2 = i1[1], i1[2], i2[1], i2[2]
            strand1, strand2 = i1[-2], i2[-2]
            diff_bp = len([_ for _, _2 in zip(seq1, seq2)
                           if _ != _2 and _ != '-' and _2 != '-'])
            text = '\t'.join([str(_) for _ in [
                             block, g1, s1, e1, strand1, g2, s2, e2, strand2, diff_bp, len(seq1)]]) + '\n'
            f_tab.write(text)
        with open(join(odir, f"{block}.aln"), 'w') as f1:
            for genome, start, end, strand, seq in info:
                f1.write(f">{genome}:{start}-{end}:{strand}\n{seq}\n")
    f_tab.close()


def seq_compare(seq_array, start, end, window_size=1000):
    # start,end is a genomic coordiante. may be reversed since some are aligned after reverse_complement
    # thus, it need to be converted prior to slicing the seqs
    assert seq_array.shape[1] == abs(end-start)
    diff_sites = np.apply_along_axis(lambda x: len(set(x)) != 1, 0, seq_array)

    if start < end:
        coord_dict = dict(zip(np.arange(0, seq_array.shape[1]),
                              np.arange(start, end, 1)))
    else:
        coord_dict = dict(zip(np.arange(0, seq_array.shape[1]),
                              np.arange(start, end, -1)))
    pos2SNP_num = {}
    for s_pos in np.arange(0, seq_array.shape[1], window_size):
        pos = s_pos+window_size/2
        _s = s_pos
        _e = s_pos+window_size
        sub_seq = diff_sites[_s:_e]
        if sub_seq.shape[0] == 0:
            SNP_num = 0
        else:
            SNP_num = sub_seq.sum()
        if pos not in coord_dict:
            # might be less than 500bp in the last section
            continue
        pos2SNP_num[coord_dict[int(pos)]] = SNP_num
    return pos2SNP_num


def get_density(block_tab, block_in_dir, ref, ofile, genome2contig2size,min_num=None):
    # reference free
    block_df, end_col, start_col, genome_col = get_fixed_block_df(
        block_tab, ref,genome2contig2size)
    block_df = block_df.set_index('Block')
    text = "ref_chr\tmid pos\tSNP\tblock_num\n"
    for block_f in tqdm(glob(f"{block_in_dir}/*.aln")):
        block = block_f.split('/')[-1].replace('.aln', '')
        records = list(SeqIO.parse(block_f, 'fasta'))
        if min_num is not None and len(records) < min_num:
            continue
        seq_list = [list(str(_.seq)) for _ in records]
        seq_list = np.array(seq_list)

        ref_block = [_
                     for _ in records
                     if _.id.split(':')[0].split('.')[0] == ref]
        ref_idx = [idx
                   for idx, _ in enumerate(records)
                   if _.id.split(':')[0].split('.')[0] == ref]
        if not ref_block:
            continue
        idx = ref_idx[0]
        seq_list = seq_list[:, seq_list[idx, :] != '-']
        contig, start, end = block_df.loc[block,  [genome_col, start_col, end_col]]
        contig = contig.split('.')[-1]
        start, end = int(start), int(end)
        pos2snp = seq_compare(seq_list, start, end, 1000)
        for pos, v in pos2snp.items():
            text += f"{contig}\t{pos}\t{v}\t{block}\n"
    with open(ofile, 'w') as f1:
        f1.write(text)




def fix_block_df(block_df,genome2contig2size):
    # since the strand would significant affect the coordinates transformation
    # For example, if the strand is -1, the start and end are the ones following reverser complemented sequence.
    # we need to fix it in advance
    block_df.loc[block_df['strand']==1,['corrected_Start','corrected_End']] = block_df.loc[block_df['strand']==1,['Start','End']].values
    block_df.loc[block_df['strand.1']==1,['corrected_Start.1','corrected_End.1']] = block_df.loc[block_df['strand.1']==1,['Start.1','End.1']].values

    def recoordinated(start,end,size):
        forward_idx = np.arange(size+1)
        reverse_idx = forward_idx[::-1]
        newstart,newend = reverse_idx[start],reverse_idx[end]
        return newstart,newend
        
    idx_l = []
    idx1_l = []
    for idx,row in tqdm(block_df.loc[block_df.isna().any(1),:].iterrows()):
        if row['strand'] == -1:
            name = row['Genome1']
            genome,c = name.split('.')
            size = genome2contig2size[genome][c]
            start,end = row['Start'],row['End']
            s,e = recoordinated(start,end,size)
            idx_l.append((idx,s,e))
        if row['strand.1'] == -1:
            name = row['Genome2']
            genome,c = name.split('.')
            size = genome2contig2size[genome][c]
            start,end = row['Start.1'],row['End.1']
            s,e = recoordinated(start,end,size)
            idx1_l.append((idx,s,e))
    block_df.loc[[_[0] for _ in idx_l],['corrected_Start','corrected_End']] = [_[1:] for _ in idx_l]
    block_df.loc[[_[0] for _ in idx1_l],['corrected_Start.1','corrected_End.1']] = [_[1:] for _ in idx1_l]
    return block_df

# 坐标转换的痛苦流程 (aln_coordinate_convertor)
# alignment block的strand为-1的，其start,end其实是完整的contig reverse complemen后的坐标，所以其实无意义，需要将其转换成完整contig上的坐标，需要用以上  fix_block_df，会得到start < end的corrected的列 (都是0-coordinate)
# 对于vcf中的SNP pos (1-coordinate)，由于其是基于concat_aln，而这个存在gap，所以需要将其转换成 无gap时的坐标 (1-coordinate)，用aln_coordinate_convertor
# 将 无gap时的坐标，map到contig上的坐标时，需要有几步
# 1. 从concat的block对应到每个block自己与基因组的位置，即找到block_df中的row
# 2. map 无gap时的坐标 (1-coordinate) 到基因组的位置
def aln_coordinate_convertor(aln, ref, pos_list):
    records = list(SeqIO.parse(aln, 'fasta'))
    seq_list = [list(str(_.seq)) for _ in records]
    # seq_array = np.array(seq_list)
    num_sites = len(seq_list[0])
    # seq_array.apply
    # get index of the reference
    ref_idx = [idx for idx, _ in enumerate(records) if _.id == ref]
    ref_idx = ref_idx[0]
    set_pos_list = set(pos_list)
    # pos before removing gap which is the same as the vcf
    
    S1 = np.arange(num_sites) + 1
    # nogap_seq_array = seq_array[:, seq_array[ref_idx, :] != '-']
    # S2 = np.arange(seq_list.shape[1]) + 1  # pos after removing
    
    ori_idx = 1
    pos2pos = []
    for aln_idx, bp in tqdm(zip(S1, seq_list[ref_idx]),
                        total=num_sites):
        # if bp =='-' and idx+1 in pos_list:
        #     print(idx+1)
        if bp == '-':
            pos2pos.append((aln_idx,'NA'))
            continue
        if aln_idx in set_pos_list:
            pos2pos.append((aln_idx,ori_idx))
        ori_idx += 1
    return pos2pos


def get_fixed_block_df(block_tab,ref,genome2contig2size):
    block_df = pd.read_csv(block_tab, sep='\t')
    block_df.loc[:, 'g1'] = [_.split('.')[0] for _ in block_df['Genome1']]
    block_df.loc[:, 'g2'] = [_.split('.')[0] for _ in block_df['Genome2']]
    block_df = fix_block_df(block_df,genome2contig2size)
    
    _block_info = {}
    for idx,row in tqdm(block_df.iterrows(),total=block_df.shape[0]):
        b = row['Block']
        if ref not in list(row[['g1','g2']]): continue
        s,e = (row['corrected_Start'],row['corrected_End']) if row['g1'] == ref else (row['corrected_Start.1'],row['corrected_End.1']) 
        if (b,s,e) not in _block_info:
            _block_info[(b,s,e)] = row
    # assert len(_block_info)==block_orders
    block_df = pd.concat(_block_info.values(),axis=1).T
    assert len(block_df['g2'].unique()) == 1
    if ref in list(block_df['g1']):
        block_df = block_df.drop(columns=['Genome2','g2']+[_ for _ in block_df.columns if _.endswith('.1')])
    else:
        block_df = block_df.drop(columns=['Genome1','g1']+[_.replace('.1','') for _ in block_df.columns if _.endswith('.1')])
    block_df.columns = [_.replace('.1','') for _ in block_df.columns]
    _m = {'Genome2':'Genome','Genome1':'Genome'}
    block_df.columns = [_m.get(_,_) for _ in block_df.columns]
    return block_df

def map_back_pos(pos2pos, block_order, block_tab, ref,genome2contig2size, ofile=None):
    block_df = get_fixed_block_df(block_tab,ref,genome2contig2size)

    block_orders = [_ for _ in open(block_order).read().split('\n') if _]
    block_df = block_df.set_index('Block')
    block_df = block_df.loc[block_orders]

    block_df.loc[:, 'real_size'] = (block_df['corrected_End'] - block_df['corrected_Start']).abs()
    block_df.loc[:, 'cumsum_pos'] = block_df['real_size'].cumsum()

    with open(ofile, 'w') as f1:
        new_SNP_coordinates = "concat pos\tref_chr\tpos\tblock\n"
        f1.write(new_SNP_coordinates)
        for concat_pos, ref_pos in tqdm(pos2pos):
            pos = ref_pos-1  # since the ref_pos is 1-coordinates
            row = block_df.loc[block_df.loc[:,'cumsum_pos'] > pos, :].iloc[0, :]
            contig = row['Genome'].split('.')[-1]
            # distance to relative to the end of the block
            dis = row['cumsum_pos'] - pos
            if row['corrected_Start'] > row['corrected_End']:  # it is reversed
                n_pos = row['corrected_End'] + dis
            else:
                n_pos = row['corrected_End'] - dis
            f1.write(f"{concat_pos}\t{contig}\t{n_pos+1}\t{row.name}\n")
            # n_pos is 1-coordinated


def summary_SNP_pos(infile, block_tab, ofile,genome2contig2size, window_size=1000):
    # start,end is a genomic coordiante
    # thus, it need to be converted prior to slicing the seqs
    snp_df = pd.read_csv(infile, sep='\t')
    ref = snp_df['ref_chr'][0].split('_')[0]
    snp_df.loc[:, 'chr'] = [int(_.split('_')[-1]) for _ in snp_df['ref_chr']]
    snp_df = snp_df.sort_values(['chr', 'pos'], ascending=True)

    block_df,end_col, start_col, genome_col = get_fixed_block_df(block_tab,ref,genome2contig2size)
    f1 = open(ofile, 'w')
    f1.write('ref_chr\tmid pos\tSNP\tblock_num\n')
    for align_chr in tqdm(block_df[genome_col].unique()):
        chr = align_chr.split('.')[-1]
        sub_df = block_df.loc[block_df[genome_col] == align_chr, :]
        for idx, row in sub_df.groupby(['Start', 'End']).head(1).iterrows():
            start, end = int(row[start_col]), int(row[end_col])
            if start <=end:
                iter_range = np.arange(start, end+1, window_size)
            else:
                iter_range = np.arange(start, end+1, window_size*-1)
            for s_pos in iter_range:
                pos = s_pos+window_size/2
                sub_snp_df = snp_df.loc[(snp_df['ref_chr'] == chr) & (
                    snp_df['pos'] < s_pos+window_size) & (snp_df['pos'] >= s_pos), :]
                f1.write(
                    f"{chr}\t{pos}\t{sub_snp_df.shape[0]}\t{row['Block']}\n")
    f1.close()
