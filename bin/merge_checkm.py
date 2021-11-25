import os
from glob import glob
from os.path import *

import click
import pandas as pd

from dating_workflow.step_script import process_path
from tqdm import tqdm

remaining_columns = [
    "marker lineage",
    "# markers",
    "# marker sets",
    "Completeness",
    "Contamination",
]


def main(indir, ofile, clean):
    indir = process_path(indir)
    ofile = process_path(ofile)

    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    all_tsv = glob(join(indir, "*", "storage", "bin_stats_ext.tsv"))
    result = {}
    for each_f in tqdm(all_tsv):
        text = open(each_f).read().strip("\n")
        gid, v = text.split("\t")
        v = eval(v)
        v = {k: _v for k, _v in v.items() if not k.startswith("GCN")}
        result[gid] = v
    new_df = pd.DataFrame.from_dict(result, orient="index")
    if clean:
        new_df = new_df.loc[:, remaining_columns]
    new_df.to_csv(ofile, sep="\t")


def sing(infile):
    result = {}
    for row in tqdm(open(infile)):
        text = row.strip("\n")
        gid, v = text.split("\t")
        v = eval(v)
        v = {k: _v for k, _v in v.items() if not k.startswith("GCN")}
        result[gid] = v


@click.command()
@click.option("-i", "indir", help="actually is the odir of checkm process")
@click.option("-o", "ofile", help="output tab file")
@click.option(
    "-clean",
    "clean",
    is_flag=True,
    default=False,
    help="clean some columns which might be meaningless when you are using the protein sequence as input",
)
def cli(indir, ofile, clean):
    main(indir, ofile, clean)


if __name__ == "__main__":
    cli()
