import os
from glob import glob
from os.path import *

import click
from Bio import SeqIO

from dating_workflow.step_script import convert_genome_ID, process_path


def main(indir, suffix, num_genomes):
    all_genes = glob(join(indir, f'*.{suffix}'))
    all_genes_name = [basename(_).replace(f".{suffix}", "")
                      for _ in all_genes]

    genes_need = []
    for f, name in zip(all_genes, all_genes_name):
        records = SeqIO.parse(f, format='fasta')
        num_records = len(list(records))
        if num_records >= num_genomes:
            genes_need.append(name)
    return genes_need


@click.command()
@click.option('-i', 'indir')
@click.option("-o", "ofile", default=None, help="path of outfile.")
@click.option('-s', 'suffix', default='aln')
@click.option("-gl", "genome_list", default=None, help="it will read 'selected_genomes.txt', please prepare the file, or indicate the alternative name or path.")
@click.option('-num_p', 'num_percentage', default=None, type=float, help="1-100")
@click.option('-num', 'num_genomes', type=int, default=None)
def cli(indir, ofile, suffix, num_percentage, num_genomes, genome_list):
    if genome_list is None:
        genome_list = join(indir, 'selected_genomes.txt')
    with open(genome_list, 'r') as f1:
        gids = f1.read().split('\n')
    gids = [convert_genome_ID(_) for _ in gids]

    if num_genomes is None:
        num_genomes = len(gids)
    if num_percentage is None:
        num_percentage = 100

    num_genomes = int(num_genomes) * int(num_percentage) / 100

    genes = main(indir, suffix, num_genomes)

    ofile = process_path(ofile)
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    print(f"found {len(genes)} meet requirement.")
    with open(ofile, 'w') as f1:
        f1.write('\n'.join(genes))


if __name__ == "__main__":
    cli()
