"""
The steps for resorting the OG_df
1. each value within the table must be a single locus instead of comma separated multiple locus
2. extract the rows which doesn't have locus on the backbone
3. reinsert extracted row from step2 into the raw_df(after removing step2 rows) based on
    a. use the leftmost locus/genome as standard to move whole row
    b. find the nearest locus (include up and down / left and right)
        c. considerate the position of current locus and its neighbours
        d. deal with flipped situation or normal order situation
    e. return the index for insertion. (it might insert into the upstream or downstream, it depends)
    f. insert it, the inserted list would gradually growth

"""

import os
import pickle
from os.path import *

import click
import numpy as np
import pandas as pd
from tqdm import tqdm

tmp_dir = join(os.environ.get('PWD'), '.tmp')


def preprocess_locus_name(locus):
    if '|' in locus:
        locus = str(locus).split('|')[-1].split(' ')[0]
    else:
        locus = str(locus).split(' ')[0]
    locus = locus.strip()
    return locus


# order_func = lambda x: int(preprocess_locus_name(x).split('_')[-1].strip(string.ascii_letters)) if not pd.isna(x) else pd.np.inf


def order_a_column(Acolumn, ordered_locus):
    Acolumn = Acolumn.dropna()
    locus2OG = {v: k for k, v in Acolumn.items()}

    sorted_column = sorted(Acolumn, key=lambda x: ordered_locus[preprocess_locus_name(x)])
    sorted_series = pd.Series(sorted_column,
                              index=[locus2OG[locus] for locus in sorted_column])
    return sorted_series


def get_next_locus_idx(genome_ordered_col, used_locus, ordered_locus):
    reorder_col = sorted(genome_ordered_col + [used_locus],
                         key=lambda x: ordered_locus[preprocess_locus_name(x)] if not pd.isna(x) else np.inf)
    # this is a ascending list
    # add this locus into these OG, and reorder it.
    # maybe reverse alignment (flipped), so it should compare the next and up idx.
    used_locus_idx = reorder_col.index(used_locus)
    if used_locus_idx == len(reorder_col) - 1:
        # the last one, no next, it mean it is the largest one
        # if right order, up_locus_idx is the last one, this locus should insert into the end
        # if reverse order, up_locus_idx is the first one, this locus should insert into the front
        next_locus_idx = -1
        # right order mean up_locus_idx is 
    else:
        next_locus = reorder_col[used_locus_idx + 1]
        next_locus_idx = genome_ordered_col.index(next_locus)

    if used_locus_idx == 0:
        # the first one, no up, it mean it is the smallest one
        # if right order, next_locus_idx is the first one, this locus should insert into the front
        # if reverse order, next_locus_idx is the last one, this locus should insert into the end
        up_locus_idx = -1
    else:
        up_locus = reorder_col[used_locus_idx - 1]
        up_locus_idx = genome_ordered_col.index(up_locus)

    if (next_locus_idx == -1 and next_locus_idx == 0) or \
            (up_locus_idx == -1 and next_locus_idx == len(genome_ordered_col)):
        # insert into the end
        return len(genome_ordered_col)
    elif (next_locus_idx == -1 and next_locus_idx == len(genome_ordered_col)) or \
            (up_locus_idx == -1 and next_locus_idx == 0):
        # insert into the front
        return 0

    if next_locus_idx < up_locus_idx:
        # if original next one smaller than up one, it mean reversed.
        # so we should insert this locus into the front of up one instead of next one
        return up_locus_idx
    else:
        return next_locus_idx


def main(infile, backbone_column_idx=0, subset_columns=None):
    OG_df = pd.read_csv(infile, sep='\t', index_col=0, low_memory=False)
    sub_idx = OG_df.index[OG_df.applymap(lambda x: ',' in str(x)).any(1)]
    # tmp
    # gfiles = glob('/mnt/home-backup/jjtao/protein/*.faa')
    # genome2order_tuple = {}
    # for g in tqdm(gfiles):
    #     gname = basename(g).replace('.faa','')
    #     if gname in OG_df.columns:
    #         genome2order_tuple[gname] = [preprocess_locus_name(_.description)
    #                                      for _ in SeqIO.parse(g, format='fasta')]

    # cal the contig cover
    # total_contig_for_each = {g:len(set([_.split(' ')[1].split(':')[0]
    #                                     for _ in gnames] ))
    #                          for g, gnames in genome2order_tuple.items()}
    # backbone_column_ori = OG_df.iloc[:, backbone_column_idx]
    # ref_OGs = OG_df.index[~backbone_column_ori.isna()]
    # for other_g in total_contig_for_each:
    #     gnames = genome2order_tuple[other_g]
    #     ref_gnames = OG_df.loc[ref_OGs,other_g]
    #     ref_gnames = [_ for _ in list(ref_gnames) if not pd.isna(_)]
    #     its_g = 'pass'

    if len(sub_idx) != 0:
        raise Exception("It contains duplicated genes within single OG. Please use `split_out_duplicated.py` first. ")
    if isinstance(backbone_column_idx, int):
        backbone_column_ori = OG_df.iloc[:, backbone_column_idx]
        used_genome = OG_df.columns[backbone_column_idx]
    else:
        backbone_column_ori = OG_df.loc[:, backbone_column_idx]
        used_genome = backbone_column_idx
    gap_OGs = OG_df.index[backbone_column_ori.isna()]
    gap_OG_df = OG_df.loc[gap_OGs, :]
    gap_OG_df = gap_OG_df.loc[~gap_OG_df.isna().all(1), :]
    genome2order_tuple = pickle.load(open(join(tmp_dir, 'genome2order_tuple'), 'rb'))
    ordered_locus = genome2order_tuple[used_genome]
    ordered_locus = [_ for v in ordered_locus for _ in v]
    # This is mainly for expanding locus located at multiple contigs into a linear ordered list
    ordered_locus = {preprocess_locus_name(locus): _
                     for _, locus in enumerate(ordered_locus)}
    backbone_c = order_a_column(backbone_column_ori, ordered_locus)
    order_OG_without_gap = list(backbone_c.index)
    order_OG_df = OG_df.reindex(backbone_c.index)
    # only backbone, without inserted genes
    tqdm.write("%s OG need to be reinserted into an ordered table" % gap_OG_df.shape[0])

    for gap_OG, row in tqdm(gap_OG_df.iterrows(),
                            total=gap_OG_df.shape[0]):
        used_genome, used_locus = [(genome, locus)
                                   for genome, locus in row.items()
                                   if not pd.isna(locus)][0]

        ordered_locus = genome2order_tuple[used_genome]
        # it return nested list
        # each list is sets of genes within a contig
        # it need to flatten this concated list
        ordered_locus = [_ for v in ordered_locus for _ in v]
        ordered_locus = {preprocess_locus_name(locus): _
                         for _, locus in enumerate(ordered_locus)}
        # iterative to get the first not nan one to order whole row.
        genome_ordered_col = list(OG_df.reindex(order_OG_without_gap).loc[:, used_genome])
        # order this genome among all order_OG_without_gap OGs. (may nan, but for backbone is all full and order.)
        insert_idx = get_next_locus_idx(genome_ordered_col,
                                        used_locus,
                                        ordered_locus)
        # find the next locus(need to justify the order of this index, not arbitrary +1), and take the index as the insert index.
        order_OG_without_gap.insert(insert_idx, gap_OG)
    final_OG_df = pd.concat([order_OG_df, gap_OG_df]).reindex(order_OG_without_gap)

    if subset_columns is not None:
        subset_columns = open(subset_columns).read().split('\n')
        diff_c = final_OG_df.columns.differences(subset_columns)
        if len(diff_c) != 0:
            raise IOError(f"{len(diff_c)} names provided in subset_columns isn't detected at your table.")
        same_c = [_ for _ in subset_columns if _ in final_OG_df.columns]
        final_OG_df = final_OG_df.loc[:, same_c]
        final_OG_df = final_OG_df.loc[~final_OG_df.isna().all(1), :]
    return final_OG_df


@click.command(help="resort the orthlogroups table. Default with the first columns as a backbone and reinsert the others orthlogroups into the centre of their neighbours")
@click.option("-i", "infile", help="input file. The file must be file after splitting duplicated OG.")
@click.option("-o", "ofile", help='output file')
@click.option("-bc", "backbone_column", help='which columns you want to taken as backbone. default is the first one', default=0)
@click.option("-sub", "subset_columns", help="subset list of genomes in the header of the table")
def cli(infile, ofile, backbone_column, subset_columns):
    final_OG_df = main(infile, backbone_column, subset_columns)
    final_OG_df.to_csv(ofile, sep='\t', index=1)


if __name__ == '__main__':
    cli()
