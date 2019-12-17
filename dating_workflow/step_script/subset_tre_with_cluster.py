"""
Script for convert cluster file into normal selected_genomes_id.txt
"""
import os
import random
from collections import defaultdict
from os.path import dirname, exists

import click


@click.command()
@click.option("-i", "infile")
@click.option("-o", "ofile")
@click.option("-s", "seed", default=100)
def cli(infile, ofile, seed):
    random.seed(seed)

    new_genomes_list = []
    genome2cluster_id = {}
    cluster2gnames = defaultdict(list)
    for row in open(infile, 'r'):
        if not row.startswith('SequenceName'):
            gname = row.split('\t')[0]
            clusterID = row.split('\t')[-1].strip('\n')
            genome2cluster_id[gname] = clusterID
            if clusterID != '-1':
                cluster2gnames[clusterID].append(gname)
            else:
                new_genomes_list.append(gname)
    print("origin contains %s genomes " % len(genome2cluster_id))
    for _, gnames in cluster2gnames.items():
        gname = random.choice(gnames)
        new_genomes_list.append(gname)

    print("remained %s genomes " % len(new_genomes_list))
    if "/" not in ofile:
        ofile = "./" + ofile
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    with open(ofile, 'w') as f1:
        f1.write('\n'.join(new_genomes_list))


if __name__ == "__main__":
    cli()
