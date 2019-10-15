import pandas as pd
from os.path import *
from Bio import SeqIO
import os
from subprocess import check_call
import plotly.express as px
homolog_dict = {'mmoA':['amoA',
                        'bmoA',
                        'pmoA',
                        'pxmA'],
                'mmoB':['amoB','bmoB','pmoB','pxmB'],
                'mmoC':['amoC','bmoC','pmoC','pxmC'],
                'OCC':['hao','haoA','nrfA','hdh','hzo','ONR','OTR'],
                'NXR_NAR_A':['nxrA','narG','dmsA','torA','unknown'],
                'NXR_NAR_B':['nxrB','narH','dmsB','torC'],
                }
in_dir = '/home-user/thliao/project/nitrogen_cycle/curated_genes'
odir = './validate'
in_file = '/home-user/thliao/project/nitrogen_cycle/manually_curated_N_cycle_genes.xlsx'

def modify_record(args):
    record,name = args
    record.id = record.id + '_' + name
    record.name = record.description = ''
    return record

from api_tools.itol_func import to_color_Clade,renamed_tree,to_node_symbol
def write2colorbranch_clade(id2info,odir,info2color,treefile, unique_id,info_name='type',**kwargs):
    content = to_color_Clade(id2info,info2color,treefile,info_name,**kwargs)
    info_name = info_name.replace('/','_')
    with open(join(odir, f'{unique_id}_{info_name}_colorbranch_clade.txt'), 'w') as f1:
        f1.write(content)
        
        
def get_id2tax(id,df):
    tmp_df = df.copy()
    tmp_df = tmp_df.set_index("AA accession")
    tmp_df.index = [str(_).strip() for _ in tmp_df.index]
    tax = tmp_df.loc[id,'phylum/class']
    return tax
    

redo = True
    
if not exists(odir):
    os.makedirs(odir)
in_df = pd.read_excel(in_file)
in_df = in_df.loc[in_df.loc[:,'note']!='removed',:]
in_df.loc[:,'phylum/class'] = in_df.loc[:,'phylum/class'].fillna('failed')
new_in_df = in_df.loc[in_df.loc[:,'phylum/class'].isin(['failed','ENV']),:]
dropped_ids = [_.strip() for _ in new_in_df.loc[:,'AA accession']]
for homolog_name,gene_list in homolog_dict.items():
    if not homolog_name.startswith('NXR_NAR'):
        continue
    odir = abspath(odir)
    total_fa = []
    for gene_name in gene_list:
        gn = gene_name
        fa_file = f'./{gn}.faa' 
        if exists(fa_file):
            records = SeqIO.parse(fa_file,format='fasta')
            records = [_ for _ in records if _.id not in dropped_ids]
            records = list(map(modify_record,
                               [(record,gn) for record in records]))
            now_ids = [_.id for _ in total_fa]
            for record in records:
                if record.id not in now_ids:
                    total_fa.append(record)
                    now_ids.append(record.id)

    with open(join(odir,f'{homolog_name}.faa'),'w') as f1:
        SeqIO.write(total_fa,f1,format='fasta-2line')
    if not exists(f'{odir}/{homolog_name}.aln') or redo:
        check_call(f'mafft --maxiterate 1000 --genafpair --thread -1 {odir}/{homolog_name}.faa > {odir}/{homolog_name}.aln', shell=1)
        
    tree_suffix = 'treefile'
    if not exists(f'{odir}/{homolog_name}.{tree_suffix}') or redo:     
        #check_call(f'FastTree {odir}/{homolog_name}.aln > {odir}/{homolog_name}.{tree_suffix}', shell=1)
        check_call(f'iqtree -nt 30 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre {odir}/{homolog_name} -s {odir}/{homolog_name}.aln', shell=1)
        
    id2info = {id:id.rpartition('_')[-1] for id in now_ids}
    colors = px.colors.qualitative.Dark24 + px.colors.qualitative.Light24
    info2color = dict(zip( list(set(id2info.values())),
                          colors))
    renamed_tree(f'{odir}/{homolog_name}.{tree_suffix}',
                 f'{odir}/{homolog_name}_new.{tree_suffix}')
    new_text = to_node_symbol( f'{odir}/{homolog_name}_new.{tree_suffix}')
    with open(join(odir, f'{homolog_name}_node_bootstrap.txt'), 'w') as f1:
        f1.write(new_text)
    write2colorbranch_clade(id2info,odir,info2color,f'{odir}/{homolog_name}_new.{tree_suffix}',homolog_name,info_name='gene cluster')
    id2info = {id:get_id2tax(id.rpartition('_')[0],in_df) for id in now_ids}
    colors = list(px.colors.qualitative.Dark24) + list(px.colors.qualitative.Light24)
    info2color = dict(zip(sorted(list(set(id2info.values()))),
                          colors))
    write2colorbranch_clade(id2info,odir,info2color,f'{odir}/{homolog_name}_new.{tree_suffix}',homolog_name,info_name='taxonomy')