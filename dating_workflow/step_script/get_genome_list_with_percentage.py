"""
read alignment/faa/trimal output for counting the number of genes contained at each genomes. And filtering out the genomes contain lower number ones.
"""
import click
from glob import glob
from os.path import *
from dating_workflow.step_script import convert_genome_ID,convert_genome_ID_rev,process_path
from Bio import SeqIO
import os
from collections import defaultdict
from tqdm import tqdm

def main(indir,suffix,num_genes):
    all_genes = glob(join(indir,f'*.{suffix}'))
    gid2num = defaultdict(int)
    
    tqdm.write('reading all genes...')
    for f in tqdm(all_genes):
        records = SeqIO.parse(f,format='fasta')
        for r in records:
            gid = convert_genome_ID_rev(r.id)
            gid2num[gid] += 1
    genomes = {k for k,v in gid2num.items() if v >=num_genes}     
    tqdm.write(f"detect {len(genomes)} match given params...")
    return genomes
    
@click.command()
@click.option('-i','indir')
@click.option("-o", "ofile", default=None, help="path of outfile.")
@click.option('-s','suffix',default='aln')
@click.option("-t","num_total_genes",default=120)
@click.option('-num_p','num_percentage',default=None,help="1-100")
@click.option('-num','num_genes',default=None)
def cli(indir,ofile,suffix,num_percentage,num_genes,num_total_genes):
    
    if num_genes is None and num_percentage is None:
        num_percentage = 100
        num_genes = num_total_genes
    elif num_genes is None:
        num_percentage = int(num_percentage)
        num_genes = num_total_genes*num_percentage/100
    else:
        num_genes = int(num_genes)
    tqdm.write(f"Filter out genomes which only contain {num_genes} ")
    genomes = main(indir,suffix,num_genes)
    
    ofile = process_path(ofile)
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    with open(ofile,'w') as f1:
        f1.write('\n'.join(list(genomes)))
    
if __name__ == "__main__":
    cli()