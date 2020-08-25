"""
This script is mainly used to run Bayestraits v3 and summarizing its result

Multistate + MCMC and likelihood ratio test

It accept 
1. input tree (maybe use .contree, .ufboot would implemented latly)
2. input metadata
3. output directory
"""
import os
from os.path import *
from subprocess import check_call
import multiprocessing as mp
import click
from ete3 import Tree

from for_software.for_bayestraits.toolkit.construct_kit import nw2nexus, get_tags
from for_software.for_bayestraits.toolkit.get_result import get_result, summaized_r,summaized_rate

# intree = './trees/iqtree/over20p_bac120.formatted.newick'
# inmetadata = './bayesTraits_test/m2nm.txt'
# odir = './bayestraits_habitat'

bt_exe = expanduser("~/software/BayesTraitsV3.0.2-Linux/BayesTraitsV3")

def run(cmd):
    check_call(cmd, shell=True)



@click.command()
@click.option('-i', 'intree',help="input tree in newick format. (should have internal node name)")
@click.option('-im', 'inmetadata',help="input metadata which contain two columns using tab as separator. No header. First line should be the id from the input tree. The second line should be the state you want to assign. ")
@click.option('-o', 'odir',help="output directory")
@click.option('-color', 'color_dict',default="M:#0000ff;N:#D68529",help="color scheme for states. default is 'M:#0000ff;N:#D68529' . ")
@click.option("-extra_cmd","extra_cmd",default='',help="extrac command for bayestraits.")
def cli(intree, inmetadata, odir,color_dict,extra_cmd):
    main(intree, inmetadata, odir,color_dict,extra_cmd,3)
    

def main(intree, inmetadata, odir,color_dict,extra_cmd,tree_format=3,threshold=None):
    if not exists(odir):
        os.makedirs(odir)
    tree_prepared_file = join(odir, basename(intree))
    metadata_pre_file = join(odir, 'metadata.txt')

    # format tree
    t = Tree(intree, format=tree_format)
    all_gids = set(t.get_leaf_names())
    new_tree_text = nw2nexus(t)
    with open(tree_prepared_file, 'w') as f1:
        f1.write(new_tree_text)

    # check metadata
    states_collect = []
    m_text = []
    for row in open(inmetadata):
        if ' ' in row:
            row = row.replace(' ','\t')
        if row.count('\t') >=2:
            print(f"maybe something wrong, please check the output {metadata_pre_file}")
        if row.split('\t')[0] in all_gids:
            m_text.append(row.strip('\n'))
            states_collect+=list(row.strip('\n').split('\t')[1])
    states = set(states_collect)
    random_states = ''.join(list(states)[:2])
    with open(metadata_pre_file, 'w') as f1:
        f1.write('\n'.join(m_text))
 
    # run
    complex_model = ["1", "2", "PriorAll exp 10", "Stones 100 1000"]
    simple_model = ["1", "2", "PriorAll exp 10",
                    f"RestrictAll q{random_states}", 
                    "Stones 100 1000"]
    tags_list = get_tags(intree)

    os.makedirs(join(odir, 'complex_m'), exist_ok=True)
    os.makedirs(join(odir, 'simple_m'), exist_ok=True)
    with open(join(odir, 'complex_m', 'params.txt'), 'w') as f1:
        f1.write('\n'.join(complex_model + tags_list +
                           [f"LF {join(odir, 'complex_m', 'bst_complex')}",
                            extra_cmd, 
                            'run']))
    with open(join(odir, 'simple_m', 'params.txt'), 'w') as f1:
        f1.write('\n'.join(simple_model + tags_list +
                           [f"LF {join(odir, 'simple_m', 'bst_simple')}",
                            extra_cmd, 
                            'run']))
    
    cmd1 = f"{bt_exe} {tree_prepared_file} {metadata_pre_file} < {join(odir, 'complex_m', 'params.txt')} > /dev/null"
    cmd2 = f"{bt_exe} {tree_prepared_file} {metadata_pre_file} < {join(odir, 'simple_m', 'params.txt')} > /dev/null"
    
    print("start to run cmd")
    cmds = [cmd1,cmd2]
    with mp.Pool(processes=2) as tp:
        r = list(tp.imap(run, cmds))
        
    if isinstance(color_dict,str):
        cat2color = color_dict.split(';')
        cat2color = {_.split(':')[0]:_.split(':')[1] for _ in cat2color if _}
        
    text = get_result(join(odir, 'complex_m', 'bst_complex.Log.txt'),
                      cat2info=cat2color)

    with open(join(odir, 'complex_habitat_prob.itol.txt'), 'w') as f1:
        f1.write(text)

    text = summaized_r(complex_f=join(odir, 'complex_m', 'bst_complex.Stones.txt'),
                       simple_f=join(odir, 'simple_m','bst_simple.Stones.txt'),
                       key='Stone')
    text2 = summaized_rate(complex_f=join(odir, 'complex_m', 'bst_complex.Log.txt'),
                           key='Iteration')
    
    with open(join(odir, 'summaized_result.txt'), 'w') as f1:
        f1.write(text+'\n'+text2)


if __name__ == '__main__':
    cli()
