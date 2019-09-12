import pandas as pd
from Bio import SeqIO
import os
from os.path import join
from tqdm import tqdm
manually_info = '../manually_curated_N_cycle_genes.xlsx'  # manually curated genes with ko info? or not
subject_info_df = pd.read_excel(manually_info)
subject_info_df = subject_info_df.set_index('AA accession')

final_df = pd.read_csv('./concated_all.csv', sep='\t', index_col=0)
target_ = 'hao_hdh_εHao (HaoA)_nrfA'
genes = target_.split('_')
os.makedirs(target_, exist_ok=1)

r_id = list(final_df.index[final_df.loc[:, 'Gene name(N metabolism)'].isin(genes)])
r_id += list(final_df.index[final_df.loc[:, 'paralog KO name'].isin(genes)])
r_id = set(r_id)
with open(join(target_, 'protein.faa'), 'w') as f1:
    records = SeqIO.parse('./first_extract_seq.faa', format='fasta')
    reads = [_ for _ in tqdm(records) if _.id in r_id]
    remained_seq = r_id.difference(set([_.id for _ in reads]))
    records = SeqIO.parse('../new_grab/output/first_extract_seq.faa', format='fasta')
    reads += [_ for _ in tqdm(records) if _.id in remained_seq]
    SeqIO.write(reads, f1, format='fasta-2line')
    f1.flush()
    _subdf = subject_info_df.loc[subject_info_df.loc[:, 'gene name'].isin(genes), :]
    for _,row in _subdf.iterrows():
        f1.write('>%s\n' % _)
        f1.write('%s\n' % row['AA sequence(seq)'])

import seaborn as sns
colors = sns.color_palette('Set1',len(genes)).as_hex()

template = open('/home-user/thliao/db/itol_template/dataset_binary_template.txt').read()

template = template.replace("FIELD_SHAPES",
                 "FIELD_SHAPES\t"+ '\t'.join(['2']*len(genes)))
template = template.replace("FIELD_LABELS",
                 "FIELD_LABELS\t"+ '\t'.join(genes))
template = template.replace("FIELD_COLORS",
                 "FIELD_COLORS\t"+ '\t'.join(colors[:len(genes)]))
for idx,g in enumerate(genes):
    _subdf = subject_info_df.loc[subject_info_df.loc[:, 'gene name']==g, :]
    label = '\t'.join(['1' if _id==idx else '-1' for _id,v in enumerate(genes)])
    for _,row in _subdf.iterrows():
        template += '%s\t%s\n' % (_,label)
with open(join(target_,'curated_genes.txt'),'w') as f1:
    f1.write(template)

#"iqtree -nt 32 -m MFP -bb 1000 -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -pre iqtree -s protein.aln"