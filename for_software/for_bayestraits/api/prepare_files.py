"""
This script is mainly used to generate tags and tree file for bayestraits.


"""
import sys
import os
from os.path import *
sys.path.insert(0,dirname(dirname(dirname(dirname(__file__)))))

import click
from ete3 import Tree

from for_software.for_bayestraits.toolkit.construct_kit import nw2nexus, get_tags

@click.command()
@click.option('-i', 'intree',help="input tree in newick format. (should have internal node name)")
@click.option('-o', 'odir',help="output directory")
def cli(intree, odir):
    main(intree, odir,)
    

def main(intree, odir,tree_format=3,):
    if not exists(odir):
        os.makedirs(odir)
    tree_prepared_file = join(odir, basename(intree))

    # format tree
    t = Tree(intree, format=tree_format)
    new_tree_text = nw2nexus(t)
    with open(tree_prepared_file, 'w') as f1:
        f1.write(new_tree_text)

    # run
    tags_list = get_tags(intree)
    with open(join(odir,'tags.txt'), 'w') as f1:
        f1.write('\n'.join(tags_list))


if __name__ == '__main__':
    cli()
