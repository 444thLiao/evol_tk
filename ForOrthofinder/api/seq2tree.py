"""
Take multiple-fasta file to construct a tree with iqtree / FastrTree

"""
import click
from subprocess import check_call
import os
from os.path import join, dirname, abspath

mafft_p = "/usr/local/bin/mafft"
iqtree_p = "/usr/local/bin/iqtree"
fasttree_p = "/home-user/software/FastTree/FastTreeMP"


def run_cmd(cmd):
    check_call(cmd,
               shell=True)


@click.command()
@click.option("-i", "infile")
@click.option("-o", "odir", default=None)
@click.option("-n","--name",default=None)
@click.option("-redo", is_flag=True, default=False)
def main(infile, odir, name,redo):
    infile = abspath(infile)
    if name is None:
        basename = os.path.basename(infile).split('.')[0]
    else:
        basename = name
    if odir is None:
        odir = f'./{basename}'
    odir = abspath(odir)
    os.makedirs(odir, exist_ok=True)
    cmd = f"{mafft_p} --maxiterate 1000 --genafpair --thread -1 {infile} > {odir}/{basename}.aln"
    if (not os.path.exists(f"{odir}/{basename}.aln")) or redo:
        run_cmd(cmd)
    # test fasttree (for quick view)
    os.makedirs(f"{odir}/quick_view", exist_ok=True)
    cmd = f"{fasttree_p} {odir}/{basename}.aln > {odir}/quick_view/{basename}.newick"
    if (not os.path.exists(f"{odir}/quick_view/{basename}.newick")) or redo:
        run_cmd(cmd)
    # iqtree
    os.makedirs(f"{odir}/iqtree/", exist_ok=True)
    cmd = f"{iqtree_p} -s {odir}/{basename}.aln -nt 10 -m MFP -pre {odir}/iqtree/{basename} -redo\n"
    if (not os.path.exists(f"{odir}/iqtree/iqtree.treefile")) or redo:
        run_cmd(cmd)


if __name__ == '__main__':
    main()
    # python3 ~/script_api/ForOrthofinder/api/seq2tree.py -i ./OG0000151.faa