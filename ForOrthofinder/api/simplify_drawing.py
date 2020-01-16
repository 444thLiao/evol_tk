
from glob import glob
import os,io
from collections import defaultdict
import pandas as pd
from os.path import *
from tqdm import tqdm
import pickle
import re
from ete3 import Tree

Nif_list = ['K02591',
            'K02586',
            'K02588']
Nod_list = ['K14658',
            'K14659',
            'K14666',
            'K09695',
            'K09694',
            ]

Nif_list = ['K02584',
 'K02585',
 'K02586',
 'K02587',
 'K02588',
 'K02589',
 'K02590',
 'K02591',
 'K02592',
 'K02593',
 'K02594',
 'K02595',
 'K02596',
 'K02597',
 'K03737',
 'K03839',
 'K04488',
 'K15790']
Nod_list = [
 'K14657',
 'K14658',
 'K14659',
 'K14660',
 'K14661',
 'K14666',
 'K12546',
 'K18904',
 'K09694',
 'K09695',
 ]

info_df = pd.read_csv("./nodnif_annotated_df.tsv",sep='\t',index_col=0)

OG_df = pd.read_csv('./subset_og.tsv',sep='\t',index_col=0)
tmp_dir = './.tmp'
if exists(join(tmp_dir,'genome2gene_info')):
    genome2gene_info = pickle.load(open(join(tmp_dir,'genome2gene_info'), 'rb'))
    genome2order_tuple = pickle.load(open(join(tmp_dir,'genome2order_tuple'), 'rb'))


## rename gbk (do only once)

def refine_tree(tree_file,indir=None):
    text = open(tree_file).read()
    new_text = re.sub("\[[0-9]+\]",'',text)
    new_text = new_text.replace('OROOT','')
    t = Tree(new_text)
    right_ordered_names = t.get_leaf_names()
    if indir is not None and exists(indir):
        right_ordered_names = [rn
                               for rn in right_ordered_names
                               if exists(join(indir,f"{rn}.fna"))]
    return right_ordered_names