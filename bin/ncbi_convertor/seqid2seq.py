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

@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id ')
@click.option('-d', 'database', default='nuccore', help='default is nuccore. ')
@click.option('-o', 'ofile', help='output fasta file')
@click.option('-f', 'format', help='default is fasta. fasta or genbank')
@click.option('-s', 'size', help='default is 100')
def cli(infile, ofile, database,format,size):
    id_list = open(infile).read().strip().split('\n')
    convertor = NCBI_convertor(id_list,database)
    if (not exists(dirname(ofile))) and '/' in ofile:
        os.system(f"mkdir -p {dirname(ofile)}")
    if format.lower()=='fasta'    :
        seqs = convertor.get_seq(batch_size=int(size))
        with open(ofile,'w') as f1:
            SeqIO.write(seqs,f1,'fasta-2line')
    elif format.lower()=='genbank':
        seqs = convertor.get_seq(batch_size=int(size),
                                 preset='genbank')
        with open(ofile,'w') as f1:
            SeqIO.write(seqs,f1,'genbank')
if __name__ == '__main__':
    cli()
