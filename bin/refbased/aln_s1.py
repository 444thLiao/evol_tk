import os
from os.path import *
from subprocess import check_call
import click
import json
import pandas as pd
from ete3 import Tree
import numpy as np
from for_software.for_EPA.parse_jplace import parse_tree_with_edges

guppy_exe = "/home-user/thliao/download/pplacer-Linux-v1.1.alpha19/guppy "


def aln(ref_newick,ref_phy,query_faa,name):
    odir = dirname(realpath(query_faa))
    cmd = f"cd {odir}; papara -a -t {ref_newick} -s {ref_phy} -q {realpath(query_faa)} -r -n {name} -j 30"
    os.system(cmd)
    cmd2 = f"epa-ng --split {ref_phy} {odir}/papara_alignment.{name} --redo"
    check_call(cmd2,shell=1)
    check_call(f"mv reference.fasta  {odir}/; mv query.fasta {odir}/;",shell=1)
    cmd3 = f"epa-ng --ref-msa {odir}/reference.fasta --tree {ref_newick} --query {odir}/query.fasta -T 30 -w {odir}  --model 'LG' --redo"
    check_call(cmd3,shell=1)
    check_call(f'mv {odir}/epa_result.jplace {odir}/epa_{name}.jplace; mv {odir}/epa_info.log {odir}/epa_{name}.log',shell=1)

def get_length(ref_newick):
    tre = Tree(ref_newick,format=3)
    n2plength = {n.name:n.dist for n in tre.traverse()}
    all_dist = list(n2plength.values())
    #np.mean(all_dist),np.max(all_dist)
    percentile95 = np.percentile(all_dist,95)
    return n2plength,percentile95

def epa(query_faa,name,n2plength,percentile95):
    odir = dirname(realpath(query_faa))
    jplace = f"{odir}/epa_{name}.jplace"
    cmd = f"{guppy_exe} to_csv {jplace} > {odir}/jplace_{name}.csv"
    check_call(cmd,shell=1)
    cmd = f"guppy edpl {jplace} --csv -o {odir}/edpl_{name}.csv"
    check_call(cmd,shell=1)
    edpl_df = pd.read_csv(f"{odir}/edpl_{name}.csv",index_col=0,header=None)
    df = pd.read_csv(f"{odir}/jplace_{name}.csv")
    df = df.sort_values('like_weight_ratio',ascending=False).groupby('name').head(1)
    df.index = df['name']
    df = df.reindex(edpl_df.index)
    df.loc[:,'EDPL'] =  edpl_df[1]
    df = df.reindex(columns=['edge_num', 'like_weight_ratio','distal_length', 'pendant_length','EDPL'])
    obj = json.load(open(jplace))
    used_tree = obj["tree"]
    tree,node_name2edge_num = parse_tree_with_edges(used_tree)
    e2n = {v:k for k,v in node_name2edge_num.items()}
    df.loc[:,'node'] = [e2n[_] for _ in df['edge_num']]
    df.loc[:,'parental length'] =  [n2plength[_] for _ in df['node']]
    final_df = df.loc[(df['pendant_length']<=percentile95) & (df['EDPL']<=0.1),:] 
    return final_df

def epa_get_data(query_faa, ref_newick, ref_phy,name):
    aln(ref_newick,ref_phy,query_faa,name)
    n2plength,percentile95 = get_length(ref_newick)
    final_df = epa(query_faa,name,n2plength,percentile95)
    return final_df




@click.command()
@click.option('-q', 'query_faa')
@click.option('-t', 'ref_newick')
@click.option('-ref', 'ref_phy',)
@click.option('-n', 'name')
def cli(query_faa, ref_newick, ref_phy,name):
    epa_get_data(query_faa, ref_newick, ref_phy,name)


if __name__ =='__main__':
    cli()

# from ete3 import Tree

# tree = Tree(ref_newick,format=3)
# name2node = {_.name:_ for _ in tree.traverse()}
# all_names = list(name2node)

# group1 = []
# for target_n in ['I50']:
#     group1+=list([_.name for _ in name2node[target_n].traverse()])
# group1.remove('I50')
# sub_df = final_df.loc[final_df['node'].isin(group1)]

# names = list(sub_df.index)
    