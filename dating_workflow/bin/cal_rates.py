"""
This script is mainly used for batch estimating the rate 
For  amino acids sequences
nucleotide


Notions:
using LG model

"""

import click
from glob import glob
from Bio import SeqIO
from ete3 import Tree
from os.path import *
import os
from bin.multiple_sbatch import sbatch_all
from bin.seqC import aln2phy
from api_tools.for_tree.format_tree import read_tree
def read_caltree(intree):
    t = open(intree).read().strip().split('\n')
    if ';' not in t[0]:
        ref_tree = Tree(t[1],8,quoted_node_names=1)
    else:
        ref_tree = read_tree(intree)
    return ref_tree
template = """ndata = 1
seqtype = 2
runmode = 0
model = 2
aaRatefile = /home-user/thliao/script/evol_tk/dating_workflow/ctl_template/lg.dat
Mgene = 0
fix_alpha = 0
alpha = 0.5
Malpha = 0
ncatG = 5
fix_rho = 1
rho = 0.
clock = 1
getSE = 1
cleandata = 0
RateAncestor = 0
"""

def prepare_cmd(aln,intree,odir):
    assert aln.endswith('aln')
    gids = aln2phy(aln,aln.replace('.aln','.phylip'))   
    gene = aln.split('/')[-1].replace('.aln','')
    print(gene,len(gids))
    ref_tree = read_caltree(intree)
    tre = ref_tree.copy()
    try:
        tre.prune(gids)
    except Exception as e:
        print(e,aln)

    for n in tre.traverse():
        if not n.is_leaf():
            n.name = ''
    tre.name = '@4.5'
    if not exists(f"{odir}/{gene}"):
        os.makedirs(f"{odir}/{gene}")
    with open(f'{odir}/{gene}/{gene}.tre','w') as f1:
        f1.write(f"{len(gids)} 1\n")
        f1.write(tre.write(format=8,format_root_node=1).replace('NoName',''))
    ctl_file = f"{odir}/{gene}/{gene}.ctl"
    with open(ctl_file,'w') as f1:
        c = template 
        c += f"seqfile = " + aln.replace('.aln','.phylip') + '\n'
        c += f"treefile = " + f'{gene}.tre' + '\n'
        c += f"outfile = " + f'{gene}.out' + '\n'
        f1.write(c)
    return f"cd {dirname(ctl_file)}; /home-user/thliao/software/paml4.7/bin/codeml {basename(ctl_file)}"


if __name__ == '__main__':
    import sys
    args = sys.argv[1:]
    cmds = []
    
    aln = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/phy_files/P39_prot/withEuk/NP_043148.1.aln'
    intree = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/cal/Plastid_B15.newick'
    odir = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/phy_files/P39_prot/withEuk/rates'
    cmd = prepare_cmd(aln,intree,odir)
    # 

    # sbatch_all(cmds,reset_workdir=1,thread_per_tasks=1,prefix_name='rate',fixed_cluster='others')



