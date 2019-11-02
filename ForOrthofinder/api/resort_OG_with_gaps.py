import pandas as pd
import click
from tqdm import tqdm


def preprocess_locus_name(locus):
    locus = str(locus).split('|')[-1].split(' ')[0]
    locus = locus.strip()
    return locus


order_func = lambda x: int(preprocess_locus_name(x).split('_')[1]) if not pd.isna(x) else pd.np.inf


def order_a_column(Acolumn):
    Acolumn = Acolumn.dropna()
    locus2OG = {v: k for k, v in Acolumn.items()}

    sorted_column = sorted(Acolumn, key=order_func)

    sorted_series = pd.Series(sorted_column,
                              index=[locus2OG[locus] for locus in sorted_column])
    return sorted_series

def get_next_locus_idx(genome_ordered_col,used_locus):
    reorder_col = sorted(genome_ordered_col + [used_locus], key=order_func)
    # this is a ascending list
    # add this locus into these OG, and reorder it.
    # maybe reverse alignment, so it should compare the next and up idx.
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
        
    if (next_locus_idx == -1 and next_locus_idx == 0) or (up_locus_idx == -1 and next_locus_idx==len(genome_ordered_col)):
        # insert into the end
        return len(genome_ordered_col)
    elif (next_locus_idx == -1 and next_locus_idx == len(genome_ordered_col)) or (up_locus_idx == -1 and next_locus_idx==0):
        # insert into the front
        return 0

    if next_locus_idx < up_locus_idx:
        # if original next one smaller than up one, it mean reversed.
        # so we should insert this locus into the front of up one instead of next one
        return up_locus_idx
    else:
        return next_locus_idx


def main(infile, backbone_column_idx=0):
    OG_df = pd.read_csv(infile, sep='\t', index_col=0)
    sub_idx = OG_df.index[OG_df.applymap(lambda x: ',' in str(x)).any(1)]
    if len(sub_idx) != 0:
        raise Exception("It contains duplicated genes within single OG. Please use `split_out_duplicated.py` first. ")
    if isinstance(backbone_column_idx, int):
        backbone_column_ori = OG_df.iloc[:, backbone_column_idx]
    else:
        backbone_column_ori = OG_df.loc[:, backbone_column_idx]
    gap_OGs = OG_df.index[backbone_column_ori.isna()]
    gap_OG_df = OG_df.loc[gap_OGs, :]
    backbone_c = order_a_column(backbone_column_ori)
    order_OG_without_gap = list(backbone_c.index)
    order_OG_df = OG_df.reindex(backbone_c.index)
    tqdm.write("%s OG need to be reinserted into an ordered table" % len(gap_OGs))

    for gap_OG, row in tqdm(gap_OG_df.iterrows(),
                            total=gap_OG_df.shape[0]):

        used_genome, used_locus = [(genome, locus)
                                   for genome, locus in row.items()
                                   if not pd.isna(locus)][0]
        # iterative to get the first not nan one to order whole row.
        genome_ordered_col = list(OG_df.reindex(order_OG_without_gap).loc[:, used_genome])
        # order this genome among all order_OG_without_gap OGs. (may nan, but for backbone is all full and order.)
        insert_idx = get_next_locus_idx(genome_ordered_col,used_locus)
        # find the next locus(need to justify the order of this index, not arbitrary +1), and take the index as the insert index.
        order_OG_without_gap.insert(insert_idx, gap_OG)
    final_OG_df = pd.concat([order_OG_df, gap_OG_df]).reindex(order_OG_without_gap)
    return final_OG_df


@click.command(help="resort the orthlogroups table. Default with the first columns as a backbone and reinsert the others orthlogroups into the centre of their neighbours")
@click.option("-i", "infile", help="input file. The file must be file after splitting duplicated.")
@click.option("-o", "ofile", help='output file')
@click.option("-bc", "backbone_column", help='which columns you want to taken as backbone. default is the first one', default=0)
def cli(infile, ofile, backbone_column):
    final_OG_df = main(infile, backbone_column)
    final_OG_df.to_csv(ofile, sep='\t', index=1)


if __name__ == '__main__':
    cli()
