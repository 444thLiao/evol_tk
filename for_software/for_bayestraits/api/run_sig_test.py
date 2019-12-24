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

import click
from ete3 import Tree

from for_software.for_bayestraits.toolkit.construct_kit import nw2nexus, get_tags
from for_software.for_bayestraits.toolkit.get_result import get_result, summaized_r,summaized_rate

# intree = './trees/iqtree/over20p_bac120.formatted.newick'
# inmetadata = './bayesTraits_test/m2nm.txt'
# odir = './bayestraits_habitat'

bt_exe = expanduser("~/software/BayesTraitsV3.0.2-Linux/BayesTraitsV3")


@click.command()
@click.option('-i', 'intree')
@click.option('-o', 'odir')
@click.option('-states',"states",default="M;N")
def main(intree, odir,states):
    if not exists(odir):
        os.makedirs(odir)
    tree_prepared_file = join(odir, basename(intree))
    metadata_pre_file = join(odir, 'metadata.txt')

    template_mode = ["1", "2", "PriorAll exp 10", "Stones 100 1000"]
    # format tree
    t = Tree(intree, format=3)
    new_tree_text = nw2nexus(t)
    if not exists(tree_prepared_file):
        with open(tree_prepared_file, 'w') as f1:
            f1.write(new_tree_text)
            
    params = []
    sig_test_dir = join(odir,'sig_test')
    for node in t.traverse():
        if not node.is_leaf():
            node_name = node.name.split('_')[0]
            tag1,tag2 = node.get_leaf_names()[0],node.get_leaf_names()[-1]
            for state in states.split(';'):
                complex_mode = template_mode+["AddTag FNode " + ' '.join([tag1,tag2]),
                                                f"Fossil Node01 FNode {state}"]
                
                ofile = join(sig_test_dir, f'{node_name}_sigtest_{state}.txt')
                with open(ofile, 'w') as f1:
                    f1.write('\n'.join(complex_mode +
                                       [f"LF {ofile.replace('.txt','')}",
                                        'Run']))
                cmd = f"{bt_exe} {tree_prepared_file} {metadata_pre_file} < {ofile}"
                if not exists(ofile.replace('.txt', '')+'.Stones.txt'):
                    params.append(cmd)
    with mp.Pool(processes=10) as tp:
        r = list(tqdm(tp.imap(run, params), total=len(params)))
                    
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
                           [f"LF {join(odir, 'complex_m', 'bst_complex')}", 'run']))
    with open(join(odir, 'simple_m', 'params.txt'), 'w') as f1:
        f1.write('\n'.join(simple_model + tags_list +
                           [f"LF {join(odir, 'simple_m', 'bst_simple')}", 'run']))

    cmd1 = f"{bt_exe} {tree_prepared_file} {metadata_pre_file} < {join(odir, 'complex_m', 'params.txt')}"
    cmd2 = f"{bt_exe} {tree_prepared_file} {metadata_pre_file} < {join(odir, 'simple_m', 'params.txt')}"

    print("start to run cmd")
    check_call(cmd1 + ' >/dev/null', shell=True)
    check_call(cmd2 + ' >/dev/null', shell=True)
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
    main()
