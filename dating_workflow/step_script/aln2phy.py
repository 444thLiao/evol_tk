import click
from glob import glob
from os.path import *
from ForOrthofinder.api.concat_aln import generate_phy_file
from Bio import AlignIO

@click.command()
@click.option('-i','infile',help='aln file')
@click.option('-o','outfile',default=None,required=None)
@click.option("-gl", "genome_list", default=None, 
              help="it will read 'selected_genomes.txt', please prepare the file, or indicate the alternative name or path.")
def cli(infile,outfile,genome_list):
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
    for infile in infiles:
        if len(infiles) != 1:
            outfile = infile.rpartition('.')[0]+'.phy'
        aln_record = AlignIO.read(infile, format='fasta')
        generate_phy_file(outfile,[(0,0,0,aln_record)],gids)
    

if __name__ == "__main__":
    cli()