"""
Format generate tree into required format
    * rooted
    * sorted
    * and also renamed the internal node for easier annotation
"""
import sys
from os.path import dirname

sys.path.insert(0, dirname(__file__))
from bin import *
import click
# required format for binary python
from api_tools.for_tree.format_tree import renamed_tree, root_tree_with
import os


def main(in_tree, o_file, outgroup_names):
    t = root_tree_with(in_tree,
                       gene_names=outgroup_names,
                       format=0)
    renamed_tree(t,
                 outfile=o_file,
                 ascending=True)


@click.command()
@click.option('-i', 'in_newick')
@click.option('-o', 'o_newick')
@click.option('-r', 'root_name', help='multiple genes could use comma to separate them. LCA would be searched and taken as outgroup')
@click.option('-f', 'force', help='overwrite?', default=False, required=False, is_flag=True)
def cli(in_newick, o_newick, root_name, force):
    if ',' in root_name:
        root_names = [_.strip() for _ in root_name.split(',')]
    else:
        root_names = [root_name.strip()]
    if not os.path.exists(dirname(o_newick)):
        os.makedirs(o_newick)

    if os.path.exists(o_newick) and not force:
        print(o_newick, ' exists and not overwrite it.')
        return
    main(in_newick, o_newick, root_names)


if __name__ == "__main__":
    cli()
