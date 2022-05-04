"""
The main purpose of this script is to retrieve sequences using seq id.
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

def pid2GI(id_list, redo=False):
    convertor = NCBI_convertor(id_list,"protein")
    suffix = 'protein2GI'
    convertor.check_cache(suffix=suffix,redo=redo)

    convertor.get_GI()
    id2gi = convertor.GI
    return id2gi

@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id ')
@click.option('-d', 'database', default='nuccore', help='default is nuccore. ')
@click.option('-o', 'ofile', help='output fasta file')
def cli(infile, ofile, database):
    id_list = open(infile).read().strip().split('\n')
    convertor = NCBI_convertor(id_list,database)
    seqs = convertor.get_seq(batch_size=100)
    if not exists(dirname(ofile)):
        os.system(f"mkdir -p {dirname(ofile)}")
    with open(ofile,'w') as f1:
        SeqIO.write(seqs,f1,'fasta-2line')
    
if __name__ == '__main__':
    cli()
