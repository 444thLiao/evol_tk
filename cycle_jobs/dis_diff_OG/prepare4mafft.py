import json
import os
from os.path import join, exists
import pandas as pd
from subprocess import check_call
from Bio import SeqIO
import seaborn as sns
# load all necessary data
odir = 'json_dump'
with open(join(odir, 'ko2og.json'), 'r') as f1:
    ko2og = json.load(f1)
with open(join(odir, 'manually_blast_r.json'), 'r') as f1:
    manually_blast_r = json.load(f1)
with open(join(odir, 'ko2og2names.json'), 'r') as f1:
    ko2og2names = json.load(f1)

og_tsv = './genome_protein_files/OrthoFinder/Results_Sep25/Orthogroups/Orthogroups.tsv'
og_df = pd.read_csv(og_tsv, sep='\t', low_memory=False, index_col=0)
genome_info = './genome_info_full.xlsx'
g_df = pd.read_excel(genome_info, index_col=0)
# generate annoating file


def generate_id2org(og_name, ofile):
    all_seq_ids = [_.id for _ in SeqIO.parse(ofile, format='fasta')]
    row = og_df.loc[og_name, :]
    id2org = {}
    for seq_id in all_seq_ids:
        org_ids = row.index[row.fillna('').str.contains(seq_id)]
        if not list(org_ids):
            print('Error', seq_id, og_name)
            continue
        org_id = org_ids[0]
        id2org[seq_id] = org_id
    return id2org
# rename
def to_label(each_og,ofile):
    template_text = open(
        '/home-user/thliao/template_txt/labels_template.txt').read()
    id2org = generate_id2org(each_og, ofile)
    full_text = template_text[::]
    for id, org in id2org.items():
        if org in g_df.index:
            name = g_df.loc[org, 'genome name']
        else:
            name = id
        full_text += '%s,%s\n' % (id, name)
    with open(join(odir, f'label_{each_og}.txt'), 'w') as f1:
        f1.write(full_text)

def to_color_strip(each_og,ofile,info_col='type'):
    template_text = open(
        '/home-user/thliao/template_txt/dataset_color_strip_template.txt').read()
    id2org = generate_id2org(each_og, ofile)
    colors = sns.color_palette('Set1').as_hex()
    
    id2info = {}
    full_text = template_text[::]
    for id, org in id2org.items():
        if org in g_df.index:
            name = g_df.loc[org, 'genome name']
            val = g_df.loc[org,info_col]
            id2info[id] = val
        else:
            name = id
    
    set_v = set(id2info.values())
    num_v = len(set_v)
    cols = colors[:num_v]
    info2col = dict(zip(set_v,cols))
    id2col = {id:info2col[info] for id,info in id2info.items()}
    annotate_text = '\n'.join(['%s,%s\n' % (id,col) for id,col in id2col.items()])
    
    legend_title = info_col
    legend_shape = ','.join(['1'] * len(info2col))
    legend_colors = ','.join([_
                              for _ in info2col.values()])
    legend_labels = ','.join(list(info2col.keys()))
    legend_text = f"""
LEGEND_TITLE,{legend_title}
LEGEND_SHAPES,{legend_shape}
LEGEND_COLORS,{legend_colors}
LEGEND_LABELS,{legend_labels}"""
    template_text = template_text.format(legend_text=legend_text,
                     dataset_label=info_col)
    with open(join(odir, f'color_{each_og}.txt'), 'w') as f1:
        f1.write(template_text+'\n'+annotate_text)

# separated OG align
og_list = [og for og_l in ko2og.values() for og in og_l]

odir = join('./align', 'single_OG')
os.makedirs(odir, exist_ok=1)
for each_og in og_list:
    fa_file = f'./genome_protein_files/OrthoFinder/Results_Sep25/Orthogroup_Sequences/{each_og}.fa'
    ofile = join(odir, each_og+'.aln')
    if not exists(ofile):
        check_call(
            f'mafft --anysymbol --thread -1 {fa_file} > {ofile}', shell=1)
    # if not exists(ofile.replace('.aln','.treefile')):
    #     check_call(f'iqtree -nt 32 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -pre {odir}/{each_og} -s {ofile}',shell=1)
    to_label(each_og,ofile)
    to_color_strip(each_og,ofile,info_col='type')
    

odir = join('./align', 'complete_ko')
os.makedirs(odir, exist_ok=1)
for ko,og_list in ko2og.items():
    
    fa_files = [f'./genome_protein_files/OrthoFinder/Results_Sep25/Orthogroup_Sequences/{each_og}.fa' for each_og in og_list]
    new_file = join(odir, ko+'.fa')
    if not exists(new_file):
        fa_file = ' '.join(fa_files)
        check_call(
            f'cat {fa_file} > {new_file}', shell=1)
    ofile = join(odir, ko+'.aln')
    if not exists(ofile):
        check_call(
            f'mafft --anysymbol --thread -1 {fa_file} > {ofile}', shell=1)
    if not exists(ofile.replace('.aln','.treefile')):
        check_call(f'iqtree -nt 32 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -pre {odir}/{ko} -s {ofile}',shell=1)
    to_label(ko,ofile)
    to_color_strip(ko,ofile,info_col='type')