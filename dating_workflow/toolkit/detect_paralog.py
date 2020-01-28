from ete3 import Tree
from glob import glob
from os.path import *
import os
from collections import defaultdict
from Bio import SeqIO
def cov_name(name):
    return 'GCA_'+name.replace('v','.').split('_')[0]


all_c_ids = open('/home-user/thliao/data/cyano_basal/rawdata/assembly_ids.list').read().split('\n')
all_c_ids = [_ for _ in all_c_ids if _]
indir = './test_paralog_cog25/82g_fasttree/'


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

[(k,_k,_v) for k,v in f2genome2gene.items() for _k,_v in v.items() if _v in potential_wrong_g]

with open("./cog25_single/84g_aln/potential_paralog.list",'w') as f1:
    f1.write('\n'.join(potential_wrong_g))


gname2diff_genome = {}
now_fa_dir = './cog25_single/82g_aln/tmp'
for seq in glob('./cog25_single/168g_aln/tmp/*.faa'):
    records = list(SeqIO.parse(seq,format='fasta'))
    name = basename(seq).replace('.faa','')
    other_fa = join(now_fa_dir,f'{name}.faa')
    o_records = list(SeqIO.parse(other_fa,format='fasta'))

    locus_same = set([_.id for _ in records]).intersection(set(_.id for _ in o_records))
    locus_diff = set([_.id for _ in records]).difference(set(_.id for _ in o_records))
    pre_same = set([_.id.split('_')[0] for _ in records]).intersection(set(_.id.split('_')[0] for _ in o_records))

    print(name,len(locus_same),len(pre_same))

    for _ in locus_diff:
        if _.split('_')[0] in pre_same:
            print(_)