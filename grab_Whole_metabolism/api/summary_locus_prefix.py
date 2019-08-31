from glob import glob
import click
from os.path import join, dirname, basename, abspath
import os
from tqdm import tqdm
from Bio import SeqIO
import pandas as pd


@click.command(help="quickly get a summary file from pokka_o")
@click.option("-i", "indir", help='input dir, normally is the output directory of prokka.')
@click.option("-o", "outfile", help='output summary file')
def main(indir, outfile, ):
    indir = abspath(indir)
    p_files = glob(join(indir, '*', "*.faa"))
    if not p_files:
        raise Exception("No faa file in %s" % indir)
    if not os.path.exists(dirname(outfile)):
        os.makedirs(dirname(outfile), exist_ok=True)

    sname2locus = {}
    for p_file in tqdm(p_files):
        sname = basename(dirname(p_file)).replace('.faa','')
        records = SeqIO.parse(p_file, format='fasta')
        record = next(records)
        locus_tag = record.id.split('_')[0]
        sname2locus[sname] = {}
        sname2locus[sname]['locus_prefix'] = locus_tag
    result_df = pd.DataFrame.from_dict(sname2locus, orient='index')
    result_df.to_csv(outfile,index_label='sample_name')
