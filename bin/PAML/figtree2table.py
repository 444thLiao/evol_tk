"""

"""
import os,sys
import re
from glob import glob
from os.path import *

from ete3 import Tree


def process_path(path):
    if not '/' in path:
        path = './' + path
    if path.startswith('~'):
        path = expanduser(path)
    if path.startswith('.'):
        path = abspath(path)
    return path

#   Tree with NHX style metadata:

#    (A:0.1,(B:0.2,(C:0.2,D:0.3):0.4[&&NHX:conf=0.01:name=NODE1]):0.5);

#     A, B, C   : leaf names
#     internal node will have the ID NODE1
#     metadata value 'conf' will be available for visualization

def sub_for(m):
    t = m.string[m.start(0):m.end(0)]
    t = t.replace(',', '-')
    t = t.replace(' ', '')
    t = t.replace('[&95%HPD={', "")
    t = t.replace('}]', "")
    return t

def get_node_name(f):
    matched_row = ''
    match = False
    for row in open(f):
        row = row.strip('\n').strip()
        if match and not row:
            break
        if row.startswith('Species tree for FigTree.  Branch lengths = posterior mean times; 95% CIs = labels'):
            match = True
            continue
        if match:
            matched_row = row
    t = Tree(matched_row.replace(' ', ''), format=8)
    for l in t.get_leaves():
        l.name = l.name.partition('_')[-1]
    return t

def main(mcmc_out_tree, out_table ):
    mcmc_out_tree_text = open(mcmc_out_tree)
    tree = None
    for row in mcmc_out_tree_text:
        if row.strip().startswith('UTREE 1 ='):
            t = row.split('UTREE 1 =')[1].strip('\n')
            t = re.sub('\[&95%HPD=.*?\]', sub_for, t)
            t = t.replace(' ', '')
            tree = Tree(t, format=1)
    if glob(join(dirname(mcmc_out_tree), '*.out')):
        outfile = glob(join(dirname(mcmc_out_tree), '*.out'))[0]
        rename_tree = get_node_name(outfile)
    else:
        rename_tree = None

    rows = ["\t".join(["node name","Posterior mean time","CIs"])]

    count = len(tree.get_leaf_names()) + 1
    for n in tree.traverse():
        if not n.is_leaf():
            dates = n.name
            if rename_tree is None:
                n.name = 'I%s' % count
            else:
                n.name = 't_n%s' % rename_tree.get_common_ancestor(n.get_leaf_names()).name
            n.add_features(ages=dates, )
            count += 1

            rows.append("\t".join([n.name,
                                   str(n.get_distance(tree.get_leaves()[0])),
                                   str(dates)]))
    if not exists(dirname(process_path(out_table))):
        os.makedirs(dirname(process_path(out_table)))
    with open(out_table,'w') as f1:
        f1.write('\n'.join(rows))


if __name__ == "__main__":
    params = sys.argv
    if not len(params) == 3:
        raise IOError("need to input path of FigTree.tre")

    main(process_path(params[1]),
         process_path(params[2]),)