#!/usr/bin/env python
"""
Format fasta-like sequences

"""
import sys
import os
from os.path import exists,join,basename,dirname
from collections import defaultdict
from glob import glob
from api_tools import to_binary_shape

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

@cli.command(help="convert multiple faa/aln/trimal (fasta-like) files to the annotations of presence and absence in itol format.")
@click.option('-i', 'indir')
@click.option('-o', 'ofile')
@click.option('-s', 'suffix',default='ffn')
@click.option('-c', 'color',default='#436bee')
def seq2itol(indir,ofile,suffix,color):
    genome2genes = defaultdict(list)
    for ffn in glob(f'{indir}/*.{suffix}'):
        gene = ffn.split('/')[-1].replace(f'.{suffix}','')
        for r in SeqIO.parse(ffn,'fasta'):
            genome2genes[r.id].append(gene)    
    text = to_binary_shape(genome2genes,same_color=color,dataset_name='Gene presence',unfilled_other=True,
                           other_params={'margin':"50"})
    with open(ofile,'w') as f1:
        f1.write(text)
if __name__ == "__main__":
    cli()
