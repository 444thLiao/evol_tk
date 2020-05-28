"""
This script is mainly for 

"""
import os
from os.path import exists, join, dirname, realpath
from subprocess import check_call

import click
import pandas as pd
from tqdm import tqdm


def run(cmd, logf):
    check_call(cmd, shell=1, stdout=logf)


def main(infile, odir, force=False, redo=False):
    if not exists(odir):
        os.makedirs(odir)
    current_dir = dirname(realpath(__file__))
    # realpath for soft link mistake.

    ext_cmd = ''
    if force:
        ext_cmd += ' -f '
    if redo:
        ext_cmd += ' -redo '

    cmd = """python3 {script} -i "{infile}" -o "{ofile}" """
    tqdm.write("executing each convetor one by one")
    with open(join(odir, 'ncbi_convertor.log'), 'w') as logf:
        script = join(current_dir, 'pid2tax.py')
        ofile = join(odir, 'pro2taxonomy.tab')
        run(cmd.format(script=script, infile=infile, ofile=ofile) + ext_cmd, logf=logf)
        script = join(current_dir, 'pid2genome.py')
        ofile = join(odir, 'pro2genome.tab')
        run(cmd.format(script=script, infile=infile, ofile=ofile) + ext_cmd, logf=logf)
        script = join(current_dir, 'pid2bio.py')
        ofile = join(odir, 'pro2Bioinfo.tab')
        run(cmd.format(script=script, infile=infile, ofile=ofile) + ext_cmd, logf=logf)
    tqdm.write("merging them. ")
    ofile = 'pro2basic_info_ALL.tab'
    df1 = pd.read_csv(join(odir, 'pro2taxonomy.tab'), sep='\t', index_col=0)
    df2 = pd.read_csv(join(odir, 'pro2genome.tab'), sep='\t', index_col=0)
    df1 = df1.reindex(df2.index)
    basic_df = pd.concat([df1, df2], axis=1)
    basic_df.to_csv(join(odir, ofile), sep='\t', index=1)
    basic_df = basic_df.loc[~basic_df.assembly_ID.isna(), :]
    ofile = 'pro2basic_info_WITH_genome.tab'
    basic_df.to_csv(join(odir, ofile), sep='\t', index=1)
    ofile = 'pro2full_info.tab'
    df3 = pd.read_csv(join(odir, 'pro2Bioinfo.tab'), sep='\t', index_col=0)
    df3.index = [_.replace("GCA_","_").replace("GCF_","_") for _ in df3.index]
    basic_df_aID = [_.replace("GCA_","_").replace("GCF_","_") for _ in basic_df.assembly_ID]
    df3 = df3.reindex(basic_df_aID)
    df3.index = basic_df.index
    full_df = pd.concat([df3, basic_df], axis=1)
    full_df = full_df.loc[:, ~full_df.isna().all(0)]
    full_df.to_csv(join(odir, ofile), sep='\t', index=1)


@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id ')
@click.option('-o', 'odir', help='output directory')
@click.option('-f', 'force', help='force overwrite?', default=False, required=False, is_flag=True)
@click.option('-redo', 'redo', help='use cache or not? default is use the cache.', default=False, required=False, is_flag=True)
def cli(infile, odir, force, redo):
    main(infile, odir, force, redo)


if __name__ == "__main__":
    cli()
