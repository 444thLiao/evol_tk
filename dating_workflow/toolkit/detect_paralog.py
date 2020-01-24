from ete3 import Tree
from glob import glob
from os.path import *
import os
from collections import defaultdict
def cov_name(name):
    return 'GCA_'+name.replace('v','.').split('_')[0]


all_c_ids = open('/home-user/thliao/data/cyano_basal/rawdata/assembly_ids.list').read().split('\n')
all_c_ids = [_ for _ in all_c_ids if _]
indir = './test_paralog_cog25/84g_fasttree/'


f2genome2gene = defaultdict(dict)
for newick in glob(join(indir,'*.newick')):
    fname = basename(newick).replace('.newick','')
    new_t = join(dirname(newick),'renamed',basename(newick))
    if not exists(dirname(new_t)):
        os.makedirs(dirname(new_t))
    t = Tree(newick)
    for l in t.get_leaves():
        f2genome2gene[fname][cov_name(l.name)] = l.name
        l.name = cov_name(l.name)
    t.write(outfile=new_t)

for newick in glob(join(indir,'renamed','*.newick')):
    t = Tree(newick)
    gname = basename(newick).replace('.newick','')
    num_non_cyano = len(set(all_c_ids).difference(set(t.get_leaf_names())))
    print(num_non_cyano,gname)


potential_wrong_g = []
for newick in glob(join(indir,'*.newick')):
    t = Tree(newick)

    leaf_names = list(t.get_leaf_names())
    for idx,n in enumerate(leaf_names[1:-1]):
        idx +=1
        genome_n = cov_name(n)
        if genome_n in all_c_ids:
            if cov_name(leaf_names[idx-1]) not in all_c_ids and cov_name(leaf_names[idx+1]) not in all_c_ids:
                potential_wrong_g.append(n)


for seq in glob(join('./seq_'))
with open