"""
Script for convert cluster file into normal selected_genomes_id.txt
"""
import os
import random
from collections import defaultdict
from os.path import dirname, exists
from ete3 import Tree
import click
from numpy.core.fromnumeric import sort

def get_most_basal_one(tree,leaves):
    """
    sort the leaves based on the ascending order
    because the traverse is following to "levelorder" strategy
    """
    names = [n.name for n in tree.traverse()]
    sorted_leaves = sorted(leaves,key=lambda x:names.index(x))
    return sorted_leaves
    

@click.command()
@click.option("-i", "infile")
@click.option("-t", "intree")
@click.option("-o", "ofile")
@click.option("-tgl", "target_genomes_list")
@click.option("-s", "seed", default=100)
def cli(infile,intree,target_genomes_list, ofile, seed):
    # random.seed(seed)
    tree = Tree(intree)
    target_genomes_set = set(target_genomes_list)
    
    uncluster_genomes = []
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
                uncluster_genomes.append(gname)
    
    cluster2repr_genome = defaultdict(list)                
    for cluster_id,leaves in cluster2gnames.items():
        sorted_leaves = get_most_basal_one(tree,leaves)
        inter = set(sorted_leaves).intersection(target_genomes_set)
        if inter:
            cluster2repr_genome[cluster_id] = list(inter)
        cluster2repr_genome[cluster_id].append(sorted_leaves[0])
        
    print("origin contains %s genomes " % len(genome2cluster_id))
    for _, gnames in cluster2repr_genome.items():
        uncluster_genomes.extend(gnames)
    new_genomes_list = list(set(uncluster_genomes))
    print("remained %s genomes " % len(new_genomes_list))
    
    
    if "/" not in ofile:
        ofile = "./" + ofile
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    with open(ofile, 'w') as f1:
        f1.write('\n'.join(new_genomes_list))


if __name__ == "__main__":
    cli()
