"""
from fasta to get id,name, description into a list
"""
import sys
from os.path import dirname
sys.path.insert(0, dirname(__file__))
from bin import *
import click
# required format for binary python
import os
from os.path import *
from Bio import SeqIO
def main(infa,ofile):
    records = SeqIO.parse(infa,format='fasta')
    if '/' not in ofile:
        ofile = './' + ofile
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    with open(ofile,'w') as f1:
        for record in tqdm(records):
            f1.write('\t'.join([record.id,record.name,record.description])+'\n')



@click.command()
@click.option('-i','infa')
@click.option('-o','ofile')
def cli(infa,ofile):
    main(infa,ofile)
        

if __name__ == "__main__":
    cli()