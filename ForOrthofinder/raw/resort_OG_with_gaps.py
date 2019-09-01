import pandas as pd
import click
from tqdm import tqdm


def preprocess_locus_name(locus):
    locus = str(locus).split('|')[-1]
    return locus


order_func = lambda x: int(preprocess_locus_name(x).split('_')[1]) if not pd.isna(x) else pd.np.inf


def order_a_column(Acolumn):
    Acolumn = Acolumn.dropna()
    locus2OG = {v: k for k, v in Acolumn.items()}

    sorted_column = sorted(Acolumn, key=order_func)

    sorted_series = pd.Series(sorted_column,
                              index=[locus2OG[locus] for locus in sorted_column])
    return sorted_series


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
    tqdm.write("%s OG need to insert into ordered table" % len(gap_OGs))

    for gap_OG, row in tqdm(gap_OG_df.iterrows(),
                            total=gap_OG_df.shape[0]):
        used_genome, used_locus = [(genome, locus)
                                   for genome, locus in row.items()
                                   if not pd.isna(locus)][0]
        genome_ordered_col = list(order_OG_df.reindex(order_OG_without_gap).loc[:, used_genome])
        reorder_col = sorted(genome_ordered_col + [used_locus], key=order_func)
        next_locus = reorder_col[reorder_col.index(used_locus) + 1]
        if not pd.isna(next_locus):
            insert_idx = genome_ordered_col.index(next_locus)
        else:
            insert_idx = len(order_OG_without_gap)
        order_OG_without_gap.insert(insert_idx, gap_OG)
    final_OG_df = pd.concat([order_OG_df, gap_OG_df]).reindex(order_OG_without_gap)
    return final_OG_df


@click.command()
@click.option("-i", "infile", help="input file. must be file after splitting duplicated.")
@click.option("-o", "ofile", help='output file')
@click.option("-bc", "backbone_column", help='which columns you want to taken as backbone. default is the first one', default=0)
def cli(infile, ofile, backbone_column):
    final_OG_df = main(infile, backbone_column)
    final_OG_df.to_csv(ofile, sep='\t', index=1)


if __name__ == '__main__':
    cli()
