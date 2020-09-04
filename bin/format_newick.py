#!/usr/bin/env python
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
from api_tools.for_tree.format_tree import renamed_tree, root_tree_with, add_cal_api, Tree, read_tree, sort_tree, earse_name, draw_cal_itol
import os
from api_tools.itol_func import to_node_symbol
from os.path import exists,join,basename
from dating_workflow.step_script import process_path

'''
ete(v3.0) toolkits format table 

      ======  ==============================================
      FORMAT  DESCRIPTION
      ======  ==============================================
      0        flexible with support values
      1        flexible with internal node names
      2        all branches + leaf names + internal supports
      3        all branches + all names
      4        leaf branches + leaf names
      5        internal and leaf branches + leaf names
      6        internal branches + leaf names
      7        leaf branches + all names
      8        all names
      9        leaf names
      100      topology only
      ======  ==============================================

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
    t = renamed_tree(in_tree,
                     outfile=o_file)
    t = root_tree_with(t,
                       gene_names=outgroup_names,
                       format=0)
    t = sort_tree(t, ascending=True)
    t.write(outfile=o_file, format=3)


def process_IO(infile, out):
    if out is None:
        out = infile.rpartition('.')[0] + '.newick'
    else:
        out = process_path(out)
        if not exists(dirname(out)):
            os.makedirs(dirname(out))
    return out


@click.group()
def cli():
    pass


@cli.command(help="simply remove all internal names of the tree")
@click.option('-i', 'in_newick')
@click.option('-o', 'out_newick', default=None)
@click.option('-f', 'tree_format', default=0)
def erase(in_newick, out_newick, tree_format):
    out_newick = process_IO(in_newick, out_newick)
    t = earse_name(in_newick, format=tree_format)
    t.write(outfile=out_newick, format=tree_format)


@cli.command(help="sort the tree. default nodes with less leaves placed bottom. ")
@click.option('-i', 'in_newick')
@click.option('-o', 'out_newick', default=None)
@click.option('-f', 'tree_format', default=0)
@click.option('-descend', 'descend', default=True, required=False, is_flag=True)
def sort(in_newick, out_newick, tree_format, descend):
    out_newick = process_IO(in_newick, out_newick)
    t = sort_tree(in_newick, ascending=descend, format=tree_format)
    t.write(outfile=out_newick, format=tree_format)


@cli.command(help="")
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


@cli.command(help="mid point root the given tree")
@click.option('-i', 'in_newick')
@click.option('-o', 'out_newick', default=None)
@click.option('-f', 'tree_format', default=0)
def mproot(in_newick, out_newick, tree_format):
    out_newick = process_IO(in_newick, out_newick)
    t = Tree(in_newick,tree_format)
    mp_n = t.get_midpoint_outgroup()
    t.set_outgroup(mp_n)
    t.write(outfile=out_newick, format=tree_format)

@cli.command()
@click.option('-i', 'in_newick')
@click.option('-o', 'out_newick', default=None)
@click.option('-f', 'tree_format', default=0)
def rename(in_newick, out_newick, tree_format):
    out_newick = process_IO(in_newick, out_newick)
    t = renamed_tree(in_newick, format=tree_format)
    with open(out_newick,'w') as f1:
        f1.write(t.write(format=3))
    # t.write(outfile=out_newick, format=new_format)


@cli.command()
@click.option('-i', 'in_newick')
@click.option('-c', 'calibration_txt')
@click.option('-o', 'out_newick', default=None)
@click.option('-f', 'tree_format', default=0)
def add_cal(in_newick, calibration_txt, out_newick, tree_format):
    out_newick = process_IO(in_newick, out_newick)
    add_cal_api(in_newick,
                out_newick=out_newick,
                calibration_txt=calibration_txt,
                format=tree_format)


@cli.command()
@click.option('-c', 'calibration_txt')
@click.option('-o', 'odir', default='./')
def itol_cal(calibration_txt, odir):
    if not exists(odir):
        os.makedirs(odir)
    draw_cal_itol(calibration_txt, odir)

@cli.command()
@click.option('-i', 'tree_file')
@click.option('-o', 'odir', default='./')
def itol_bp(tree_file, odir):
    if not exists(odir):
        os.makedirs(odir)
    text = to_node_symbol(tree_file)
    with open(join(odir,basename(tree_file).rpartition('.')[0]+'.bp.txt'),'w') as f1:
        f1.write(text)


@cli.command()
@click.option('-i', 'in_newick')
@click.option('-i2', 'in_newick2')
@click.option('-o', 'out_newick')
@click.option('-f', 'tree_format', default=0)
@click.option('-f_to', 'new_format', default=0)
@click.option('-r', 'replace', default=None, help='the leaf should be placed at i2(latter one) instead of i')
def cat(in_newick, in_newick2, out_newick, tree_format, new_format, replace):
    if replace is None:
        t = Tree()
        t1 = read_tree(in_newick, format=tree_format)
        t2 = read_tree(in_newick2, format=tree_format)
        t.add_child(t1)
        t.add_child(t2)
    else:
        t = read_tree(in_newick2, format=tree_format)
        t1 = read_tree(in_newick, format=tree_format)
        t.remove_child([_ for _ in t.children if _.name == replace][0])
        t.add_child(t1)
    text = t.write(format=new_format)


    with open(out_newick, 'w') as f1:
        f1.write(text.replace('NoName', '').replace("):0",")I0_S100:0"))


@cli.command()
@click.argument('infiles', nargs=-1)
@click.option('-o', 'out_newick')
@click.option('-f', 'tree_format', default=0)
@click.option('-f_to', 'new_format', default=0)
def mcat(infiles, out_newick, tree_format, new_format):
    """
    format_newick.py mcat a b c d -o ./test.newick
    """
    final_t = Tree()
    tree_boxs = [final_t]
    total_len = len(infiles)
    for _, f in enumerate(infiles):
        try:
            t = read_tree(f, format=tree_format)
        except IOError:
            t = Tree(f'{f};')
        tree_boxs[-1].add_child(t)
        if _ + 2 != total_len:
            # left more than one
            _t = Tree()
            if _ + 1 == total_len:
                # last one
                continue
            tree_boxs[-1].add_child(_t)
            tree_boxs.append(_t)
    text = final_t.write(format=new_format)
    with open(out_newick, 'w') as f1:
        f1.write(text.replace('NoName', ''))


@cli.command()
@click.option('-i', 'in_newick')
@click.option('-o', 'out_newick', default=None)
@click.option('-f', 'tree_format', default=0)
@click.option('-f_to', 'new_format', default=0)
def reformat(in_newick, out_newick, tree_format, new_format):
    out_newick = process_IO(in_newick, out_newick)
    t1 = read_tree(in_newick, format=tree_format)
    text = t1.write(out_newick, format=new_format)
    with open(out_newick, 'w') as f1:
        f1.write(text.replace('NoName', ''))

# preset for post analysis for the contree generate by iqtree
# 1. renamed it
# 2. reroot
# 3. sort it
@cli.command()
@click.option('-i', 'in_newick')
@click.option('-o', 'out_newick', default=None)
@click.option('-f', 'tree_format', default=0)
@click.option('-f_to', 'new_format', default=3)
@click.option('-r', 'root_name', help='multiple genes could use comma to separate them. LCA would be searched and taken as outgroup')
@click.option('-descend', 'descend', default=True, required=False, is_flag=True)
def set1(in_newick, out_newick, tree_format, new_format,root_name,descend):
    out_newick = process_IO(in_newick, out_newick)
    t = renamed_tree(in_newick, format=tree_format)
    if ',' in root_name:
        root_names = [_.strip() for _ in root_name.split(',')]
    else:
        root_names = [root_name.strip()]
    t = root_tree_with(t, gene_names=root_names)
    t = sort_tree(t, ascending=descend)
    t.write(outfile=out_newick, format=new_format)

if __name__ == "__main__":
    cli()
