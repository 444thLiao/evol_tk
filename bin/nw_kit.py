"""
Format generate tree into required format
    * rooted
    * sorted
    * and also renamed the internal node for easier annotation
"""
import sys
from os.path import dirname

sys.path.insert(0, dirname(__file__))
import click
from api_tools.for_tree.format_tree import *
import os
from os.path import exists,abspath,expanduser

_t = '''
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


def process_path(path):
    if path is None:
        return path
    if not '/' in path:
        path = './' + path
    path = abspath(expanduser(path))
    return path

def process_IO(infile, out):
    infile = process_path(infile)
    if out is not None:
        out = infile.rpartition('.')[0] + '.newick'
    elif isinstance(out,str):
        if not exists(dirname(out)):
            os.makedirs(dirname(out))
    elif out is None:
        return out
    else:
        raise IOError
    out = process_path(out)
    return out


@click.group(invoke_without_command=True)
@click.option('-f', 'print_tree_format', default=False, required=False, is_flag=True)
def cli(print_tree_format):
    if print_tree_format:
        print(_t)


@cli.command()
@click.option('-i', 'in_newick')
@click.option('-o', 'out_newick', default=None)
@click.option('-f', 'tree_format', default=0)
def erase(in_newick, out_newick, tree_format):
    out_newick = process_IO(in_newick, out_newick)
    t = earse_name(in_newick, format=tree_format)
    t.write(outfile=out_newick, format=tree_format)


@cli.command()
@click.option('-i', 'in_newick')
@click.option('-o', 'out_newick', default=None)
@click.option('-f', 'tree_format', default=0)
@click.option('-descend', 'descend', default=True, required=False, is_flag=True)
def sort(in_newick, out_newick, tree_format, descend):
    out_newick = process_IO(in_newick, out_newick)
    t = sort_tree(in_newick, ascending=descend, format=tree_format)
    t.write(outfile=out_newick, format=tree_format)


@cli.command()
@click.option('-i', 'in_newick')
@click.option('-o', 'out_newick', default=None)
@click.option('-f', 'tree_format', default=0)
@click.option('-r', 'root_name', help='multiple genes could use comma to separate them. LCA would be searched and taken as outgroup')
def reroot(in_newick, out_newick, tree_format, root_name):
    out_newick = process_IO(in_newick, out_newick)
    if ',' in root_name:
        root_names = [_.strip() for _ in root_name.split(',')]
    else:
        root_names = [root_name.strip()]
    t = root_tree_with(in_newick, gene_names=root_names, format=tree_format)
    t.write(outfile=out_newick, format=tree_format)


@cli.command()
@click.option('-i', 'in_newick')
@click.option('-o', 'out_newick', default=None)
@click.option('-f', 'tree_format', default=0)
def rename(in_newick, out_newick, tree_format):
    out_newick = process_IO(in_newick, out_newick)
    t = renamed_tree(in_newick, format=tree_format)
    t.write(outfile=out_newick, format=tree_format)


@cli.command()
@click.option('-i', 'in_newick')
@click.option('-c', 'calibration_txt')
@click.option('-o', 'out_newick', default=None)
@click.option('-f', 'tree_format', default=0)
def add_cal(in_newick, calibration_txt, out_newick, tree_format):
    out_newick = process_IO(in_newick, out_newick)
    text = add_cal(in_newick,
                   out_newick=out_newick,
                   calibration_txt=calibration_txt,
                   format=tree_format)
    with open(out_newick, 'w') as f1:
        f1.write(text)


@cli.command()
@click.option('-i', 'in_newick')
@click.option('-i2', 'in_newick2')
@click.option('-o', 'out_newick')
@click.option('-f', 'tree_format', default=0)
def cat(in_newick, in_newick2, out_newick, tree_format):
    t = Tree()
    t1 = read_tree(process_path(in_newick), format=tree_format)
    t2 = read_tree(process_path(in_newick2), format=tree_format)
    t.add_child(t1)
    t.add_child(t2)
    t.write(out_newick, format=tree_format)


@cli.command()
@click.option('-i', 'in_newick')
@click.option('-o', 'out_newick', default=None)
@click.option('-f', 'tree_format', default=0)
@click.option('-f_to', 'new_format', default=0)
def reformat(in_newick, out_newick, tree_format, new_format):
    out_newick = process_IO(in_newick, out_newick)
    t1 = read_tree(tree_format, format=tree_format)
    t1.write(outfile=out_newick, format=new_format)


if __name__ == "__main__":
    cli()
