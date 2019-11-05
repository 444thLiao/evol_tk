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
under_nodes = {'nxrA':['I49_S100', 'I137_S88', 'I88_S99', 'I71_S100'],
              'hao':['I49_S100','I32_S100'],
              'amoA':['I43_S100','I30_S100',
                      #'I15_S76'
                      ],
              'amoB':['I44_S89','I71_S100',
                      #'I14_S100','I16_S100'
                      ],
              'amoC':['I12_S100','I57_S59'],
              }

all_files = ['nr_retrieve_amoB/with_genome_Bacteria_drop_NC10_intact.faa_aln.dir/iqtree.treefile',
             'nr_retrieve_amoC/with_genome_Bacteria_drop_NC10_intact.faa_aln.dir/iqtree.treefile',
             'with_genome_amoA/with_genome_Bacteria_drop_NC10_intact.faa_aln.dir/iqtree.treefile',
             'nr_retrieve_hao/with_genome_Bacteria_intact.faa_aln.dir/iqtree.treefile',
             'nr_retrieve_nxrA/with_genome_Bacteria_drop_NC10_intact_lt_600.faa_aln.dir/iqtree.treefile'
             ]
# check_call('python3 ~/script/evolution_relative/global_search/reannotate_tree.py '+ ' '.join(all_files),shell=1)
from collections import defaultdict
g_genomes = defaultdict(list)
genome2id = []
reformat_id2ori = {}
manuall_class_ids = []
for f in all_files:
    g = f.split('/')[0].split('_')[-1]
    tree = glob(join(f,'*.sorted.newick'))[0]
    t = Tree(tree,format=1)
    all_ids = [reformat(_).strip() for _ in t.get_leaf_names()]
    if under_nodes.get(g,[]):
        _ns = [_ 
               for _ in t.iter_descendants() 
               if _.name in under_nodes.get(g,[])]
        for _n in _ns:
            manuall_class_ids += list(_n.get_leaf_names())
        manuall_class_ids = [reformat(_).strip() for _ in manuall_class_ids]
    reformat_id2ori.update(dict(zip(all_ids,
                                    [_ for _ in t.get_leaf_names()])))
    pro2genome = pd.read_csv(join(f,'info_dir','pro2genome.tab'),
                             sep='\t',index_col=0)
    pro2genome = pro2genome.fillna('None')
    for gene in all_ids:
        if gene in pro2genome.index:
            _cache = pro2genome.loc[gene,'assembly_ID']
            if isinstance(_cache,str) and _cache !='None':
                genome2id.append((_cache,gene))
                if gene in manuall_class_ids:
                    g_genomes[g].append(_cache)
            elif str(_cache) == 'None':
                print(gene,g)
            else:
                if gene in manuall_class_ids:
                    g_genomes[g] += list(_cache)
                genome2id += [(_g,gene) for _g in list(_cache)]
        else:
            print(gene,g)

genome2genes = defaultdict(list)
for g,genomes in g_genomes.items():
    for genome in genomes:
        if genome !='None':
            genome2genes[genome].append(g)

id2genes = {}
for genome,genes in genome2genes.items():
    ids = [id 
           for g,id in genome2id 
           if g == genome]
    id2genes.update({id:list(set(genes)) 
                     for id in ids})

info2style = {'amoA':{'color':'#ff0000',
                      'info':'amoA'},
              'amoB':{'color':'#ff0000',
                      'info':'amoB'},
              'amoC':{'color':'#ff0000',
                      'info':'amoC'},
              'hao':{'color':'#b68100',
                      'info':'hao'},
              'nxrA':{'color':'#4b85c1',
                      'info':'nxrA'},
              }
#odir = ''
for f in all_files:
    g = f.split('/')[0].split('_')[-1]
    tree = glob(join(f,'*.sorted.newick'))[0]
    t = Tree(tree,format=1)
    all_ids = [reformat(_).strip() for _ in t.get_leaf_names()]
        
    this_id2genes = {reformat_id2ori[_]:id2genes.get(_,[]) 
                     for _ in all_ids 
                     if _ in id2genes}
    all_text = to_binary_shape(this_id2genes,
                    info2style,
                    info_name='genes set',manual_v=['amoA','amoB','amoC','hao','nxrA'])
    with open(join(f,'gene_set.txt'),'w') as f1:
        f1.write(all_text)
