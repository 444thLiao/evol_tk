import pandas as pd
from os.path import *
import os
from tqdm import tqdm
from subprocess import check_call
import multiprocessing as mp
from glob import glob


gene_presence_tab = "./protein_annotations/kegg_diamond.crosstab"
basic_habitat_txt = "./bayestraits_habitat/over20p_bac120_fasttree/metadata.txt"
intree = "./bayestraits_habitat/over20p_bac120_fasttree/over20p_bac120.formatted.newick"
odir = './bayestraits_habitat/over20p_bac120_fasttree/gene_test'

habitat_text = open(basic_habitat_txt).read()
all_gids = [_.split('\t')[0] for _ in habitat_text.split('\n')]


gid2habitat = dict([_.split('\t') for _ in habitat_text.split('\n')])
genes_df = pd.read_csv(gene_presence_tab,sep='\t',index_col=0)


habitat_mapping_dict = {"M":'1',
                        "N":'0',
                        "NM":'-'}


for ko in tqdm(genes_df.columns):
    gid2gene = genes_df.loc[:,ko].to_dict()
    gid2gene = {k:1 if _ is not None else 0 for k,_ in gid2gene.items()}
    g_dir = join(odir,'each_gene',ko.split(':')[-1])
    if not exists(g_dir):
        os.makedirs(g_dir)
    metadata_txt = []
    for gid,v in gid2habitat.items():
        hv = habitat_mapping_dict[v]
        gene_v = str(int(gid2gene[gid]))
        metadata_txt.append(f"{gid}\t{hv}\t{gene_v}")
    with open(join(g_dir,'metadata.txt'),'w') as f1:
        f1.write('\n'.join(metadata_txt))
        
        