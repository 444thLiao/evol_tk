import os
from glob import glob
from os.path import expanduser, basename
from subprocess import check_call

from ete3 import NCBITaxa
from tqdm import tqdm

from api_tools.itol_func import *

ncbi = NCBITaxa()

blastp_pth = '/home-user/software/blast/latest/bin/blastp'


def run(cmd):
    check_call(cmd, shell=1)


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
os.makedirs(tmp_dir, exist_ok=True)
collected_gs = expanduser('~/project/nitrogen_cycle/curated_genes/')
genome2collect_genes = defaultdict(list)
gene2locus_id = defaultdict(list)

genes = ['hzsA', 'hzsB', 'hzsC', 'hao', 'hdh', 'hzo', ]
genes = ["amoA","amoB","amoC","cycA","cycB","hao","nxrA"]

do_genes = glob(join(collected_gs, '*.faa'))
do_genes = [_ for _ in do_genes if basename(_).replace('.faa', '') in genes]
for fa in tqdm(do_genes):
    gene_name = basename(fa).replace('.faa', '').strip()
    otab = join(tmp_dir, basename(fa).replace('.faa', '').strip() + '.tab')
    if (not exists(otab)) or redo:
        cmd = f"{blastp_pth} -query '{fa}' -db {db_faa} -out {otab} -num_threads 20 -outfmt 6 -evalue 1e-3"
        run(cmd)
    for row in open(otab):
        if float(row.strip('\n').split("\t")[-2]) <=1e-20:
            genome_name = row.split("\t")[1].split('_')[0]
            locus_id = row.split("\t")[1]
            genome2collect_genes[genome_name].append(gene_name)
            gene2locus_id[gene_name].append(locus_id)
genome2collect_genes = {k: list(set(v)) for k, v in genome2collect_genes.items()}


info2style = {'amoA': {'color': '#ff0000',
                       'info': 'amoA'},
              'amoB': {'color': '#ff0000',
                       'info': 'amoB'},
              'amoC': {'color': '#ff0000',
                       'info': 'amoC'},
              'hao': {'color': '#b68100',
                      'info': 'hao'},
              'cycA': {'color': '#b68100',
                       'info': 'cycA'},
              'cycB': {'color': '#b68100',
                       'info': 'cycB'},
              'nxrA': {'color': '#4b85c1',
                       'info': 'nxrA'},
              }

# for nitrification
requested_genes = ['amoA', 'amoB', 'amoC', 'hao', 'cycA', 'cycB', 'nxrA']
g2id = {f"GCA_{k.replace('v', '.')}": v for k, v in genome2collect_genes.items()}
g2id = {k: v for k, v in g2id.items() if set(v).intersection(set(requested_genes))}
from api_tools.itol_func import *
text = to_binary_shape(g2id,info2style,manual_v=requested_genes)
with open("./itol_txt/gene_set.txt","w") as f1:
    f1.write(text)

# annotate taxonomy


## special for annotate nxrA
nxrA_ids = open('/home-user/thliao/project/nitrogen_cycle/fetch_genes/query_result/manual_assigned_gids/nxrA_manual_assigned_ids').read().split('\n')
for g, ids in g2id.items():
    if 'nxrA' in ids and g not in nxrA_ids:
        ids.remove('nxrA')
g2id = {k: v for k, v in g2id.items() if set(v).intersection(set(requested_genes))}
all_text = to_binary_shape(g2id,
                           info2style,
                           info_name='genes set',
                           manual_v=requested_genes)
with open(join('./itol_txt', 'gene_set.txt'), 'w') as f1:
    f1.write(all_text)
