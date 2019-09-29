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

og_tsv = './genome_protein_files/OrthoFinder/Results_Sep27/Orthogroups/Orthogroups.tsv'
og_df = pd.read_csv(og_tsv, sep='\t', low_memory=False, index_col=0)
def rename(x):
    if pd.isna(x):
        return x
    if ', ' not in x:
        return str(x).split(' ')[0]
    else:
        return ', '.join([str(_).split(' ')[0] for _ in x.split(', ')])
og_df = og_df.applymap(rename)

genome_info = './genome_info_full.xlsx'
g_df = pd.read_excel(genome_info, index_col=0)

# special_annotate_file (out group )

# generate annoating file
def generate_id2org(og_names, ofile):
    all_seq_ids = [_.id for _ in SeqIO.parse(ofile, format='fasta')]
    if isinstance(og_names,str):
        og_names = [og_names]
    id2org = {}
    for og_name in og_names:
        row = og_df.loc[og_name, :]
        for seq_id in all_seq_ids:
            org_ids = row.index[row.fillna('').str.contains(seq_id)]
            if not list(org_ids):
                #print('Error', seq_id, og_name)
                continue
            org_id = org_ids[0]
            id2org[seq_id] = org_id
    return id2org
# rename
def to_label(each_og,ofile,ko_name=None):
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
    if ko_name is not None:
        each_og = ko_name
    with open(join(odir, f'label_{each_og}.txt'), 'w') as f1:
        f1.write(full_text)

color_scheme = {'type':{'NOB': '#e41a1c', 'comammox': '#edc31d', 
                        'AOB': '#bad5b9', 'AOA': '#358f0f'},
                'phylum/class':{'Thaumarchaeota': '#358f0f',
                                'Nitrospirae': '#edc31d',
                                'Gammaproteobacteria': '#78fce0',
                                'chloroflexi': '#e41a1c',
                                'Betaproteobacteria': '#956cb4',
                                'Alphaproteobacteria': '#8c613c'}

                }



def to_color_strip(each_og,ofile,info_col='type',ko_name=None):
    template_text = open(
        '/home-user/thliao/template_txt/dataset_color_strip_template.txt').read()
    id2org = generate_id2org(each_og, ofile)
    colors = sns.color_palette('Set1').as_hex()
    id2info = {}
    for id, org in id2org.items():
        if org in g_df.index:
            #name = g_df.loc[org, 'genome name']
            val = g_df.loc[org,info_col]
            id2info[id] = val

    set_v = set(id2info.values())
    num_v = len(set_v)
    cols = colors[:num_v]
    if info_col in color_scheme:
        info2col = color_scheme[info_col]
    else:
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
    if ko_name is not None:
        each_og = ko_name
    info_col = info_col.replace('/','_')
    with open(join(odir, f'color_{each_og}_{info_col}.txt'), 'w') as f1:
        f1.write(template_text+'\n'+annotate_text)

# ref or outgroup seq, additionally add to
ref_file = 'outgroup_total.xlsx'
ref_df = pd.read_excel(ref_file,index_col=None)
ref_df = ref_df.loc[ref_df.loc[:,'note']!='removed',:] 
def get_add_text(sub_df):
    t_text = '' 
    for _,row in sub_df.iterrows():
        aa_id = row['AA accession']
        gene_name = row['gene name']
        seq = row['seq']
        t_text+= f'>{aa_id}_{gene_name}\n{seq}\n'
    return t_text

def annotate_outgroup(sub_df):
    template_text = open(
        '/home-user/thliao/template_txt/dataset_binary_template.txt').read()
    template_text = template_text.format(filed_shape='3',
                                         filed_label='outgroup/reference')
    annotate_text = ''
    for _,row in sub_df.iterrows():
        aa_accession = row['AA accession']
        if row['type'] == 'outgroup':
            annotate_text += '\t'.join([aa_accession,'1']) +'\n'
        else:
            annotate_text += '\t'.join([aa_accession,'0']) +'\n'
    return template_text+annotate_text

# separated OG align
odir = join('./align', 'single_OG')
os.makedirs(odir, exist_ok=1)
for ko, og_list in ko2og.items():
    sub_ref_df = ref_df.loc[ref_df.loc[:,'outgroup for which KO']==ko,:]
    add_text = get_add_text(sub_ref_df)
    annotate_text = annotate_outgroup(sub_ref_df)
    for each_og in og_list:
        fa_file = f'./genome_protein_files/OrthoFinder/Results_Sep27/Orthogroup_Sequences/{each_og}.fa'
        new_fa_file = join(odir, each_og+'.fa')
        with open(new_fa_file,'w') as f1:
            f1.write(add_text+open(fa_file,'r').read())
        ofile = join(odir, each_og+'.aln')
        if not exists(ofile):
            check_call(
                f'mafft --anysymbol --thread -1 {new_fa_file} > {ofile}', shell=1)
        if not exists(ofile.replace('.aln','.treefile')):
            check_call(f'iqtree -nt 32 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre {odir}/{each_og} -s {ofile}',shell=1)
        to_label(each_og,ofile)
        to_color_strip(each_og,ofile,info_col='type')
        to_color_strip(each_og,ofile,info_col='phylum/class')

odir = join('./align', 'complete_ko')
os.makedirs(odir, exist_ok=1)
for ko,og_list in ko2og.items():
    sub_ref_df = ref_df.loc[ref_df.loc[:,'outgroup for which KO']==ko,:]
    add_text = get_add_text(sub_ref_df)
    annotate_text = annotate_outgroup(sub_ref_df)
    fa_files = [f'./genome_protein_files/OrthoFinder/Results_Sep27/Orthogroup_Sequences/{each_og}.fa' for each_og in og_list]
    new_file = join(odir, ko+'.fa')
    if not exists(new_file):
        fa_file = ' '.join(fa_files)
        check_call(
            f'cat {fa_file} > {new_file}', shell=1)
    ori_text = open(new_file,'r').read()
    with open(new_file,'w') as f1:
        f1.write(add_text+ori_text)
    
    ofile = join(odir, ko+'.aln')
    if not exists(ofile):
        check_call(
            f'mafft --anysymbol --thread -1 {new_file} > {ofile}', shell=1)
    if not exists(ofile.replace('.aln','.treefile')):
        check_call(f'iqtree -nt 32 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre {odir}/{ko} -s {ofile}',shell=1)
    to_label(og_list,ofile,ko_name=ko)
    to_color_strip(og_list,ofile,info_col='type',ko_name=ko)
    to_color_strip(og_list,ofile,info_col='phylum/class',ko_name=ko)
    with open(join(odir,f'marker_{ko}_outgroup_ref.txt'),'w') as f1:
        f1.write(annotate_text)