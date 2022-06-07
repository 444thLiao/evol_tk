"""
The main purpose of this script is to retrieve sequences information using seq id.
Only obtaining sequences using protein/nuccore assession is currently supported.
"""

try:
    from bin.ncbi_convertor import NCBI_convertor
except ModuleNotFoundError:
    import sys
    sys.path.insert(0,'/home-user/thliao/script/evol_tk')
    from bin.ncbi_convertor import NCBI_convertor
    
import click
from Bio import SeqIO
import os
from os.path import *
import pandas as pd

@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id ')
@click.option('-d', 'database', default='nuccore', help='default is nuccore. ')
@click.option('-o', 'ofile', help='output fasta file')
@click.option('-noseq', 'without_seq', is_flag=True, default=False,help='removing the seq')
def cli(infile, ofile, database,without_seq):
    id_list = open(infile).read().strip().split('\n')
    convertor = NCBI_convertor(id_list,database)
    
    if (not exists(dirname(ofile))) and '/' in ofile:
        os.system(f"mkdir -p {dirname(ofile)}")
    dbsummary = convertor.get_db_summary()
    if without_seq:
        for k in dbsummary:
            dbsummary[k].pop('sequence')
    for k,v in dbsummary.items():
        for _k,_v in v.items():
            v[_k] = str(_v).replace('\n',' ')
    final_df = pd.DataFrame.from_dict(dbsummary).T
    final_df.to_csv(ofile,sep='\t')
if __name__ == '__main__':
    cli()
