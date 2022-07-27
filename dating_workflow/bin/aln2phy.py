"""
This script is mainly for generating a phy file from a alignment file
You should pass a genome id list(most exactly, the sequence id list) for let the code know how many genomes it actually has.
"""
from glob import glob
from os.path import *

import click
from Bio import AlignIO,SeqIO

from dating_workflow.bin.concat_aln import convert_genome_ID_rev,generate_phy_file

def aln2phy(aln,phy):
    r = {_.id:_ for _ in SeqIO.parse(aln,'fasta')}
    num_s = len(r)
    num_sites = len(list(r.values())[0].seq)
    with open(phy,'w') as f1:
        f1.write(f"{num_s} {num_sites}\n")
        for k,v in r.items():
            f1.write(f"{k}               {v.seq}\n")
            
    
    
@click.command()
@click.option('-i', 'infile', help='aln file')
@click.option('-o', 'outfile', default=None, required=None)
@click.option("-gl", "genome_list", default=None,
              help="it will read 'selected_genomes.txt', please prepare the file, or indicate the alternative name or path.")
@click.option("-rm_I", "remove_identical", is_flag=True, default=False)
@click.option("-no_fill", "fill_gaps", is_flag=True, default=True)
@click.option("-name_convert", "name_convert", is_flag=True, default=True)
def cli(infile, outfile, genome_list,remove_identical,fill_gaps,name_convert):
    if not '/' in infile:
        infile = f'./{infile}'
    if '*' in infile:
        infiles = glob(infile)
    else:
        infiles = [infile]
        outfile = outfile
        # if len(infiles) >1, outfile will be useless.
    indir = dirname(infiles[0])
    if genome_list is None:
        genome_list = join(indir, 'selected_genomes.txt')
    with open(genome_list, 'r') as f1:
        gids = f1.read().split('\n')
        if name_convert:
            gids = [convert_genome_ID_rev(_) for _ in gids if _]
        gids = set(gids)
    for infile in infiles:
        if len(infiles) != 1:
            outfile = infile.rpartition('.')[0] + '.phy'
        aln_record = AlignIO.read(infile, format='fasta')
        generate_phy_file(outfile, 
                          [(0, 0, 0, aln_record)], 
                          gids,
                          fill_gaps=fill_gaps, 
                          remove_identical=remove_identical,
                          name_convertor=lambda x:convert_genome_ID_rev(str(x).split('_')[0]))


if __name__ == "__main__":
    cli()
