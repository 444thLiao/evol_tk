import pandas as pd
import sys
from os.path import exists, join,expanduser,basename
from tqdm import tqdm
from api_tools.itol_func import *
from glob import glob
from ete3 import Tree
import plotly.express as px
from ete3 import NCBITaxa
from global_search.classification_script import _classificated
from subprocess import check_call
import os

ncbi = NCBITaxa()

blastp_pth = '/home-user/software/blast/latest/bin/blastp'
def run(cmd):
    check_call(cmd,shell=1)

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


from collections import defaultdict
db_faa = './concat.faa'
redo = False
tmp_dir = join('./tmp/contained_genes/')
os.makedirs(tmp_dir,exist_ok=True)
collected_gs = expanduser('~/project/nitrogen_cycle/curated_genes/')
genome2collect_genes = defaultdict(list)
for fa in tqdm(glob(join(collected_gs,'*.faa'))):
    gene_name = basename(fa).replace('.faa','').strip()
    otab = join(tmp_dir,basename(fa).replace('.faa','').strip()+'.tab')
    if (not exists(otab)) or redo:
        cmd = f"{blastp_pth} -query '{fa}' -db {db_faa} -out {otab} -num_threads 20 -outfmt 6 -evalue 1e-50"
        run(cmd)
    for row in open(otab):
        genome_name = row.split("\t")[1].split('_')[0]
        genome2collect_genes[genome_name].append(gene_name)
genome2collect_genes = {k:set(v) for k,v in genome2collect_genes.items()}


info2style = {'amoA':{'color':'#ff0000',
                      'info':'amoA'},
              'amoB':{'color':'#ff0000',
                      'info':'amoB'},
              'amoC':{'color':'#ff0000',
                      'info':'amoC'},
              'hao':{'color':'#b68100',
                      'info':'hao'},
              'cycA':{'color':'#b68100',
                      'info':'cycA'},
               'cycB':{'color':'#b68100',
                      'info':'cycB'},
              'nxrA':{'color':'#4b85c1',
                      'info':'nxrA'},
              }

requested_genes = ['amoA','amoB','amoC','hao','cycA','cycB','nxrA']
g2id = {f"GCA_{k.replace('v','.')}":v for k,v in genome2collect_genes.items()}
g2id = {k:v for k,v in g2id.items() if set(v).intersection(set(requested_genes))}

## special for annotate nxrA
nxrA_ids = open('/home-user/thliao/project/nitrogen_cycle/fetch_genes/query_result/manual_assigned_gids/nxrA_manual_assigned_ids').read().split('\n')
for g,ids in g2id.items():
    if 'nxrA' in ids and g not in nxrA_ids:
        ids.remove('nxrA')
g2id = {k:v for k,v in g2id.items() if set(v).intersection(set(requested_genes))}
all_text = to_binary_shape(g2id,
                info2style,
                info_name='genes set',
                manual_v=requested_genes)
with open(join('./itol_txt','gene_set.txt'),'w') as f1:
    f1.write(all_text)
    
    