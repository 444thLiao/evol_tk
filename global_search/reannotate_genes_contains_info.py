import pandas as pd
import sys
from os.path import exists, join
from tqdm import tqdm
from api_tools.itol_func import *
from glob import glob
from ete3 import Tree
import plotly.express as px
from ete3 import NCBITaxa
from global_search.classification_script import _classificated
from subprocess import check_call
ncbi = NCBITaxa()

def reformat(s):
    a = s.split('_')[-1]
    if not '_' in s:
        return s
    try:
        float(a)
        return s
    except:
        if len(s.rpartition('_')[-1]) == 1:
            return s
        else:
            return s.rpartition('_')[0]

all_files = ['nr_retrieve_amoB/with_genome_Bacteria_drop_NC10_intact.faa_aln.dir/iqtree.treefile',
             'nr_retrieve_amoC/with_genome_Bacteria_drop_NC10_intact.faa_aln.dir/iqtree.treefile',
             'with_genome_amoA/with_genome_Bacteria_drop_NC10_intact.faa_aln.dir/iqtree.treefile',
             'nr_retrieve_hao/with_genome_Bacteria_intact.faa_aln.dir/iqtree.treefile',
             #'nr_retrieve_nxrA/with_genome_Bacteria_intact.faa_aln.dir/cluster_95_aln.dir/iqtree.treefile'
             ]
# check_call('python3 ~/script/evolution_relative/global_search/reannotate_tree.py '+ ' '.join(all_files),shell=1)
from collections import defaultdict
g_genomes = defaultdict(list)
genome2id = []
for f in all_files:
    
    g = f.split('/')[0].split('_')[-1]
    tree = glob(join(f,'*.sorted.newick'))[0]
    t = Tree(tree,format=1)
    all_ids = t.get_leaf_names()
    all_ids = [reformat(_).strip() for _ in all_ids]
    pro2genome = pd.read_csv(join(f,'info_dir','pro2genome.tab'),
                             sep='\t',index_col=0)
    pro2genome = pro2genome.fillna('None')
    genomes_assembly = []
    for _ in all_ids:
        if _ in pro2genome.index:
            _cache = pro2genome.loc[_,'assembly_ID']
            if isinstance(_cache,str) and _cache !='None':
                genomes_assembly.append(_cache)
                genome2id.append((_cache,_))
            elif str(_cache) == 'None':
                print(_,g)
            else:
                genomes_assembly+=list(_cache)
                genome2id+=[(_g,_) for _g in _cache]
        else:
            print(_,g)
    g_genomes[g] = genomes_assembly

genome2genes = defaultdict(list)
for g,genomes in g_genomes.items():
    for genome in genomes:
        if genome !='None':
            genome2genes[genome].append(g)

id2genes = {}
for genome,genes in genome2genes.items():
    ids = [id for g,id in genome2id if g == genome]
    id2genes.update({id:list(sorted(set(genes))) for id in ids})

info2style = {'amoA':{'color':'#ff0000',
                      'info':'amoA'},
              'amoB':{'color':'#ff0000',
                      'info':'amoB'},
              'amoC':{'color':'#ff0000',
                      'info':'amoC'},
              'hao':{'color':'#b68100',
                      'info':'hao'},
              }
odir = ''
for f in all_files:
    tree = glob(join(f,'*.sorted.newick'))[0]
    t = Tree(tree,format=1)
    all_ids = t.get_leaf_names()
    all_ids = [reformat(_).strip() for _ in all_ids]
    this_id2genes = {_:id2genes.get(_,[]) 
                     for _ in all_ids 
                     if _ in id2genes}
    all_text = to_binary_shape(this_id2genes,
                    info2style,info_name='genes set')
    with open(join(f,'gene_set.txt'),'w') as f1:
        f1.write(all_text)
