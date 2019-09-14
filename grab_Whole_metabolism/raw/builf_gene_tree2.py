import pandas as pd
from Bio import SeqIO
import os
from os.path import join
from tqdm import tqdm

# locus2info = '../new_grab/locus_info.csv'
# locus2info_df = pd.read_csv(locus2info,sep='\t',index_col=0)
final_df = pd.read_csv('./concated_all.csv', sep='\t', index_col=0)

manually_info = '../manually_curated_N_cycle_genes.xlsx'  # manually curated genes with ko info? or not
subject_info_df = pd.read_excel(manually_info)
subject_info_df = subject_info_df.set_index('AA accession')
#final_df = pd.read_csv('./concated_all.csv', sep='\t', index_col=0)
target_ = 'pmoB_amoB'
genes = target_.split('_')
os.makedirs(target_, exist_ok=1)

from collections import defaultdict
ids = []
g2locus = defaultdict(list)
for _ in genes:
    _eachs = list(final_df.index[final_df.loc[:,'Gene name(N metabolism)'].str.contains(_)])
    ids += _eachs
    for _each in _eachs:
        g2locus[_].append(_each)
ids = set(ids)
with open(join(target_, 'protein.faa'), 'w') as f1:
    records = SeqIO.parse('./first_extract_seq.faa', format='fasta')
    reads = [_ for _ in tqdm(records) if _.id in ids]
    remained_seq = ids.difference(set([_.id for _ in reads]))
    records = SeqIO.parse('../new_grab/output/first_extract_seq.faa', format='fasta')
    reads += [_ for _ in tqdm(records) if _.id in remained_seq]
    SeqIO.write(reads, f1, format='fasta-2line')

    f1.flush()
    _sub_df = subject_info_df.loc[subject_info_df.loc[:, 'gene name'].isin(genes), :]
    for _,row in _sub_df.iterrows():
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
    for l in g2locus[g]:
        template += '%s\t%s\n' % (l, label)

with open(join(target_,'curated_genes.txt'),'w') as f1:
    f1.write(template)
