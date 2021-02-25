"""
major for pack up the metadata of project/nitrification/reference_genome
"""

import json
import os
from os.path import join, exists,basename,dirname
from glob import glob
import pandas as pd
from subprocess import check_call
from Bio import SeqIO
import seaborn as sns
from ete3 import Tree
import plotly.express as px
import io
import multiprocessing as mp
from collections import defaultdict

odir = 'json_dump_v2'
indir = 'genome_protein_files_more'
og_tsv = f'./{indir}/OrthoFinder/Results_Oct01/Orthogroups/Orthogroups.tsv'
genome_info = './genome_info_full.xlsx'
MAG_annotate_file = '/home-user/thliao/data/metagenomes/update_0928_nitrification/confirmed_locus2info.tsv'
ref_file = '/home-user/thliao/project/nitrogen_cycle/nitrification/reference_genomes/outgroup and reference.xlsx'


# load orthofinder csv, and rename it
og_df = pd.read_csv(og_tsv, sep='\t', low_memory=False, index_col=0)
def rename(x):
    if pd.isna(x):
        return x
    if ', ' not in x:
        return str(x).split(' ')[0]
    else:
        return ', '.join([str(_).split(' ')[0] for _ in x.split(', ')])
    
og_df = og_df.applymap(rename)
g_df = pd.read_excel(genome_info, index_col=0)
dropped_g_df = g_df.loc[g_df.loc[:,'used']=='no',:]
for gid in dropped_g_df.index:
    if gid in og_df.columns:
        og_df.drop(gid,axis=1,inplace=True)
        
        
# name2dirname
# for mag 
name2dirname = {}
for ofile in glob(join(indir,'*.faa')):
    real_ofile = os.path.realpath(ofile)
    name = basename(ofile).replace('.faa','')
    if '/prokka_o/' in real_ofile:
        dir_name = basename(dirname(real_ofile)).replace('.faa','')
        name2dirname[name] = dir_name


color_scheme = {'type':{'NOB': '#e41a1c', 'comammox': '#edc31d', 
                        'AOB': '#bad5b9', 'AOA': '#358f0f'},
                'phylum/class':{'Thaumarchaeota': '#358f0f',
                                'Nitrospirae': '#edc31d',
                                'Gammaproteobacteria': '#78fce0',
                                'Chloroflexi': '#e41a1c',
                                'Betaproteobacteria': '#956cb4',
                                'Alphaproteobacteria': '#8c613c',
                                'Actinobacteria':'#11FF11',
                                'Planctomycetes':'#FF66bb',
                                }

                }

mag2habitat_file = '/home-user/thliao/project/nitrogen_cycle/nitrification/reference_genomes/MAG2habitat.xlsx'
metadata_for_17_parks = '/mnt/home-backup/thliao/metagenomes/17_parks/pack_up_metadata.tsv'
ofile = '/home-user/thliao/project/nitrogen_cycle/nitrification/reference_genomes/mag2info.csv'

MAG_annotate_df = pd.read_csv(MAG_annotate_file,sep='\t',index_col=0)
derep_df = MAG_annotate_df.drop_duplicates('sample name')
phylum_from_metadata_count = derep_df.groupby('phylum(from metadata)').count().iloc[:,0]
sname2phylum_metadata = dict(zip(MAG_annotate_df.loc[:,'sample name'],
                                 MAG_annotate_df.loc[:,'phylum(from metadata)']))
sname2class_metadata = dict(zip(MAG_annotate_df.loc[:,'sample name']
                                ,MAG_annotate_df.loc[:,'class(from metadata)']))
sname2phylum_metadata = {k:v for k,v in sname2phylum_metadata.items() if not pd.isna(v)}
sname2info = {k:v if v !='Proteobacteria' else sname2class_metadata[k] for k,v in sname2phylum_metadata.items()}
sname2info = {k:v for k,v in sname2info.items() if not pd.isna(v)}

mag2habitat_df = pd.read_excel(mag2habitat_file,index_col=0,)
mag2habitat = dict(zip(mag2habitat_df.index,
                       mag2habitat_df.loc[:,'classify as ']))
derep_df.loc[:,'habitat (manual)'] = [mag2habitat[_] 
                                    for _ in  list(derep_df.loc[:,'source project'])]

metadata_for_17_parks_df = pd.read_csv(metadata_for_17_parks,sep='\t',index_col=0)
sub_df = derep_df.loc[derep_df.loc[:,'source project']=='17_parks',:]
collect_list = []
for _,row in sub_df.iterrows():
    gnome_id = '_'.join(row['sample name'].split('_')[:2])
    habitat = metadata_for_17_parks_df.loc[gnome_id,'habitat (manual)']
    collect_list.append(habitat)
derep_df.loc[sub_df.index,'habitat (manual)'] = collect_list
sname2habitat = dict(zip(list(derep_df.loc[:,'sample name']),
                         list(derep_df.loc[:,'habitat (manual)'])))
sname2habitat = {k:'terrestrial' if v =='soil' else v
                    for k,v in sname2habitat.items()}
sname2habitat = {k:'freshwater' if v =='fresh water' else v
                    for k,v in sname2habitat.items()}
mag2infos = defaultdict(dict)
mag2infos['phylum/class'].update(sname2info)
mag2infos['habitat (manual)'].update(sname2habitat)
mag2info_df = pd.DataFrame.from_dict(mag2infos)
mag2info_df.to_csv(ofile,sep='\t',index=1,index_label='sample name')
####### 
# for reference genome
