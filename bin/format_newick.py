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
from api_tools.for_tree.format_tree import renamed_tree, root_tree_with,sort_tree
import os


'''
ete(v3.0) toolkits format table 
0	flexible with support values	((D:0.723274,F:0.567784)1.000000:0.067192,(B:0.279326,H:0.756049)1.000000:0.807788);
1	flexible with internal node names	((D:0.723274,F:0.567784)E:0.067192,(B:0.279326,H:0.756049)B:0.807788);
2	all branches + leaf names + internal supports	((D:0.723274,F:0.567784)1.000000:0.067192,(B:0.279326,H:0.756049)1.000000:0.807788);
3	all branches + all names	((D:0.723274,F:0.567784)E:0.067192,(B:0.279326,H:0.756049)B:0.807788);
4	leaf branches + leaf names	((D:0.723274,F:0.567784),(B:0.279326,H:0.756049));
5	internal and leaf branches + leaf names	((D:0.723274,F:0.567784):0.067192,(B:0.279326,H:0.756049):0.807788);
6	internal branches + leaf names	((D,F):0.067192,(B,H):0.807788);
7	leaf branches + all names	((D:0.723274,F:0.567784)E,(B:0.279326,H:0.756049)B);
8	all names	((D,F)E,(B,H)B);
9	leaf names	((D,F),(B,H));
100	topology only	((,),(,));

'''

def main(in_tree, o_file, outgroup_names):
    t= renamed_tree(in_tree,
                 outfile=o_file,
                 ascending=True)
    t= root_tree_with(t,
                       gene_names=outgroup_names,
                       format=0)
    t = sort_tree(t,ascending=True)
    t.write(outfile=o_file, format=3)


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
