import os
from glob import glob
from os.path import *

import click
import pandas as pd

from dating_workflow.step_script import process_path


def main(indir, ofile):
    indir = process_path(indir)
    ofile = process_path(ofile)

    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    all_tsv = glob(join(indir, '*', 'storage', 'bin_stats_ext.tsv'))
    result = {}
    for each_f in all_tsv:
        text = open(each_f).read().strip('\n')
        gid, v = text.split('\t')
        v = eval(v)
        v = {k: _v for k, _v in v.items() if not k.startswith('GCN')}
        result[gid] = v
    new_df = pd.DataFrame.from_dict(result, orient='index')
    new_df.to_csv(ofile, sep='\t')


@click.command()
@click.option('-i', 'indir', help='actually is the odir of checkm process')
@click.option('-o', 'ofile', help='output tab file')
def cli(indir, ofile):
    main(indir, ofile)


if __name__ == "__main__":
    cli()
