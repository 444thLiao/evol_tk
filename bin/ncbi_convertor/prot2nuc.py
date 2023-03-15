"""
The main purpose of this script is to convert protein id to nuccore ID


not finished
"""

try:
    from bin.ncbi_convertor import NCBI_convertor,edl
    from api_tools.third_party.metadata_parser import parse_elink_xml
except ModuleNotFoundError:
    import sys
    sys.path.insert(0,'/home-user/thliao/script/evol_tk')
    from bin.ncbi_convertor import NCBI_convertor,edl
    from api_tools.third_party.metadata_parser import parse_elink_xml
    
import click
from Bio import SeqIO
import os
from os.path import *

convertor = NCBI_convertor(id_list, "assembly")

@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id ')
@click.option('-d', 'database', default='nuccore', help='default is nuccore. ')
@click.option('-o', 'ofile', help='output fasta file')
@click.option('-f', 'format', help='default is fasta. fasta or genbank',default='fasta')
@click.option('-s', 'size', help='default is 100')
def cli(infile, ofile, database,format,size):

    id_list = open(infile).read().strip().split('\n')
    results, failed = edl.elink(
        dbfrom="protein",
        db="nuccore",
        ids=id_list,
        idtype='acc'    ,
        result_func=lambda x: parse_elink_xml(x),
    )
    pid2nids = dict(results)
    
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
