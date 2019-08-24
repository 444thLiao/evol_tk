from glob import glob
import click
import multiprocessing as mp
from subprocess import check_call
from os.path import join,dirname,basename,abspath
import os
from tqdm import tqdm
from Bio import SeqIO
import pandas as pd
@click.command()
@click.option("-i","indir")
@click.option("-o","outfile")
def main(indir,outfile,):
    indir = abspath(indir)
    p_files = glob(join(indir,'*',"*.faa"))
    if not p_files:
        raise Exception("No faa file in %s" % indir)
    if not os.path.exists(dirname(outfile)):
        os.makedirs(dirname(outfile),exist_ok=True)

    sname2locus = {}
    for p_file in tqdm(p_files):
        sname = basename(dirname(p_file))
        records = SeqIO.parse(p_file,format='fasta')
        record = next(records)
        locus_tag = record.id.split('_')[0]
        sname2locus[sname] = {}
        sname2locus[sname]['locus_prefix'] = locus_tag
    result_df = pd.DataFrame.from_dict(sname2locus,orient='index')
    result_df.to_csv(outfile)