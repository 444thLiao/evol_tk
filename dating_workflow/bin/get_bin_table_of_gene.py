"""
Construct a table with presence or absence from a list of fasta file ...
The name of the file would be taken as the gene name.
Such as
it could parse nirK.faa nirS.faa , etc.

"""
import pandas as pd
import click
from dating_workflow.step_script import convert_genome_ID_rev,convert_genome_ID
from Bio import SeqIO
from glob import glob
from os.path import *
from collections import defaultdict


def main(infiles):

    _df_dict = defaultdict(lambda :defaultdict(int))
    for f in infiles:
        name = basename(f).rpartition('.')[0]
        records = SeqIO.parse(f,format='fasta')
        records = list(records)
        for _ in records:
            _df_dict[name][convert_genome_ID_rev(_.id)] = 1
    df = pd.DataFrame.from_dict(_df_dict)
    df = df.fillna(0)
    df = df.reindex(columns=sorted(df.columns,key=lambda x:df[x].sum(),reverse=True))
    df = df.sort_values(list(df.columns))
    return df

@click.command()
@click.option("-i","indir",help="input directory which have sequence to be processed")
@click.option("-s","suffix",help="suffix of files need to be processed")
@click.option("-o","ofile",help="path of output file")
def cli(indir,suffix,ofile):
    infiles = glob(join(indir,f'*.{suffix}'))
    df = main(infiles)
    if ofile.endswith('xlsx'):
        df.to_excel(ofile)
    else:
        df.to_csv(ofile,sep='\t')

if __name__ == '__main__':
    cli()