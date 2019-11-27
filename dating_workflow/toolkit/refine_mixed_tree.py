from sklearn.from Bio import SeqIO
from ete3 import Tree,PhyloTree

# two function for dating workflow
def convert_genome_ID(genome_ID):
    # for GCA_900078535.2
    # it will return
    return genome_ID.split('_')[-1].replace('.', 'v')

def convert_genome_ID_rev(genome_ID):
    # for 900078535v2
    # it will return
    
    return 'GCA_' + genome_ID.replace('v', '.')


intree = './223556/223556.newick'
new_tree = intree.replace('.newick','_renamed.newick')
t = Tree(intree)
for l in t.get_leaves():
    l.name = convert_genome_ID_rev(l.name) +';'+l.name
text = t.write()
with open(new_tree,'w') as f1:
    f1.write(text)
    
import itertools
import pandas as pd
from tqdm import tqdm
from collections import defaultdict
t = Tree(intree)
# all_g = set([convert_genome_ID_rev(_.split('_')[0]) for _ in t.get_leaf_names()])
all_ids = t.get_leaf_names()

id_dict = defaultdict(dict)
for g1 in tqdm(t.get_leaves()):
    for g2 in t.get_leaves():
        id_dict[g1.name][g2.name] = t.get_distance(g1,g2)
dis = pd.DataFrame.from_dict(id_dict)

from sklearn.cluster import KMeans 
kmeans = KMeans(n_clusters=2, random_state=0,precompute_distances =True,tol=1e-10).fit(dis.values)
kmeans.labels_
    
id2info = defaultdict(list)
for idx,id in enumerate(dis.index):
    new_name = convert_genome_ID_rev(id.split('_')[0]) +'_'+id
    id2info[new_name] = [str(kmeans.labels_[idx])]
from api_tools.itol_func import *
text = to_binary_shape(id2info,{'1':{},'0':{}})
with open('../itol_txt/separate_tmp.txt','w') as f1:
    f1.write(text)


t = PhyloTree(intree)
# t.set_outgroup(t.get_midpoint_outgroup())
t.set_species_naming_function(lambda node: convert_genome_ID_rev(node.name.split('_')[0]) )

print(t.get_ascii(attributes=["name", "species"], show_internal=False ))
t2 = t.collapse_lineage_specific_expansions()

ntrees, ndups, sptrees =  t2.get_speciation_trees()
sptrees = list(sptrees)
print("Found %d species trees and %d duplication nodes" %(ntrees, ndups))
for spt in sptrees:
   print(len(spt.get_leaf_names()))
