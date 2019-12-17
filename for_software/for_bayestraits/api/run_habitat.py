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
from for_software.for_bayestraits.toolkit.get_result import get_result, summaized_r

# intree = './trees/iqtree/over20p_bac120.formatted.newick'
# inmetadata = './bayesTraits_test/m2nm.txt'
# odir = './bayestraits_habitat'

bt_exe = expanduser("~/software/BayesTraitsV3.0.2-Linux/BayesTraitsV3")


@click.command()
@click.option('-i', 'intree')
@click.option('-im', 'inmetadata')
@click.option('-o', 'odir')
def main(intree, inmetadata, odir):
    if not exists(odir):
        os.makedirs(odir)
    tree_prepared_file = join(odir, basename(intree))
    metadata_pre_file = join(odir, 'metadata.txt')

    # format tree
    t = Tree(intree, format=3)
    all_gids = set(t.get_leaf_names())
    new_tree_text = nw2nexus(t)
    with open(tree_prepared_file, 'w') as f1:
        f1.write(new_tree_text)

    # check metadata
    m_text = [_.strip('\n')
              for _ in open(inmetadata)
              if _.split('\t')[0] in all_gids]
    with open(metadata_pre_file, 'w') as f1:
        f1.write('\n'.join(m_text))

    # run
    complex_model = ["1", "2", "PriorAll exp 10", "Stones 100 1000"]
    simple_model = ["1", "2", "PriorAll exp 10",
                    "Restrict qNM qMN", "Stones 100 1000"]
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

    text = get_result(join(odir, 'complex_m', 'bst_complex.Log.txt'))

    with open(join(odir, 'complex_habitat_prob.itol.txt'), 'w') as f1:
        f1.write(text)

    text = summaized_r(complex_f=join(odir, 'complex_m', 'bst_complex.Stones.txt'),
                       simple_f=join(odir, 'simple_m',
                                     'bst_simple.Stones.txt'),
                       key='Stone')

    with open(join(odir, 'summaized_result.txt'), 'w') as f1:
        f1.write(text)


if __name__ == '__main__':
    main()
