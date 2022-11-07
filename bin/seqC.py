#!/usr/bin/env python
"""
Format fasta-like sequences

"""
import sys
import os
from os.path import exists,join,basename,dirname

sys.path.insert(0, dirname(__file__))
import click
from Bio import SeqIO


@click.group()
def cli():
    pass

@cli.command(help="disable soft-wrap")
@click.option('-i', 'infile')
@click.option('-o', 'ofile', default=None)
def m2l(infile, ofile):
    r = SeqIO.parse(infile,'fasta')
    with open(ofile,'w') as f1:
        SeqIO.write(r,f1,'fasta-2line')
        
        
@cli.command(help="convert aln to phylip file.")
@click.option('-i', 'infile')
@click.option('-o', 'ofile')
def aln2phy(infile,ofile):
    r = {_.id:_ for _ in SeqIO.parse(infile,'fasta')}
    num_s = len(r)
    num_sites = len(list(r.values())[0].seq)
    with open(ofile,'w') as f1:
        f1.write(f"{num_s} {num_sites}\n")
        for k,v in r.items():
            f1.write(f"{k}               {v.seq}\n")
            
            
if __name__ == "__main__":
    cli()
