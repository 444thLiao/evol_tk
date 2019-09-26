
from subprocess import check_call
from glob import glob
from os.path import join, dirname, exists, basename
import os
from tqdm import tqdm
import pandas as pd
from collections import defaultdict

ref_id = './ref_id'
target_file = './nitrification_list'
genome_info = './genome_info_full.xlsx'
odir = './annotated'

# validatation
all_g_ids = [_ for _ in open(ref_id, 'r').read().split('\n') if _]
g_df = pd.read_excel(genome_info, index_col=0)
assert not set(all_g_ids).difference(g_df.index)


sname2ko2locus = dict()
kid_list = [_ for _ in open(target_file, 'r').read().split('\n') if _]
for ofile in tqdm(glob(join(odir, '*.out'))):
    sname = basename(ofile).replace('.out', '')
    if sname in all_g_ids:
        locus2ko = dict([(_.split('\t')[0], _.split('\t')[1:])
                         for _ in open(ofile).read().split('\n')])
        ko2locus = dict()
        for kid in kid_list:
            locus = [l for l, ks in locus2ko.items() if kid in ks]
            ko2locus[kid] = locus
        sname2ko2locus[sname] = ko2locus


# check
def judge(ko2locus, type_g):
    if type_g == 'comammox':
        return all([bool(v) for v in ko2locus.values()])
    elif type_g == 'NOB':
        # nxrB K00371 nxrA K00370
        return all([bool(ko2locus.get(k)) for k in ['K00371', 'K00370']])
    elif type_g == 'AOB' or type_g == 'AOA':
        # amoA: K10944
        # amoB: K10945
        # amoC: K10946
        # hao: K10535
        # return all([bool(ko2locus.get(k)) for k in ['K10946','K10945','K10944','K10535']])
        return all([bool(ko2locus.get(k)) for k in ['K10944', 'K10535']])


confirmed_g = []
failed_g = defaultdict(list)
# get missed annotation and prepare to reannotate it
for g_id in set(all_g_ids):
    type_g = g_df.loc[g_id, 'type']
    name = g_df.loc[g_id, 'genome name']
    status = judge(sname2ko2locus[g_id], type_g)
    if status:
        confirmed_g.append(g_id)
    else:
        failed_g[type_g].append(g_id)


# auto implement first?
og_tsv = './genome_protein_files/OrthoFinder/Results_Sep25/Orthogroups/Orthogroups.tsv'
og_df = pd.read_csv(og_tsv, sep='\t', low_memory=False, index_col=0)
sub_df = og_df.reindex(columns=all_g_ids)
ko2og = defaultdict(list)
for g_id in confirmed_g:
    ko2locus = sname2ko2locus[g_id]
    for ko, locus_l in ko2locus.items():
        if locus_l:
            for l in locus_l:
                ogs = sub_df.index[sub_df.loc[:,
                                              g_id].str.contains(l).fillna(False)]
                if not list(ogs):
                    print(g_id, ko, locus_l)
                    break
                og = list(ogs)[0]
                if og not in ko2og[ko]:
                    ko2og[ko].append(og)
                break

revised_id = []
for key in list(failed_g.keys()):
    failed_g_ids = failed_g.pop(key)
    for failed_g_id in failed_g_ids:
        ko2locus = sname2ko2locus[failed_g_id]
        for k, locus_l in ko2locus.items():
            if not locus_l:
                backup_id = og_df.loc[ko2og[ko], failed_g_id]
                if len(list(backup_id)) != 1:
                    continue
                backup_id = list(backup_id)[0]
                if pd.isna(backup_id):
                    continue
                sname2ko2locus[failed_g_id][k] = [backup_id]
                revised_id.append((failed_g_id, k, backup_id))
        if not judge(sname2ko2locus[failed_g_id], g_df.loc[failed_g_id, 'type']):
            failed_g[key].append(failed_g_id)

odir = './reannotate'
os.makedirs(odir,exist_ok=True)
for type_g, g_ids in tqdm(failed_g.items()):
    if type_g == 'comammox':
        db_list = ['amoA','amoB','amoC','hao','haoA','HAO'] + ['nxrA','nxrB']
    elif type_g.lower() == 'aob' or type_g.lower() == 'aoa':
        db_list = ['amoA','amoB','amoC','hao','haoA','HAO']
    elif type_g.lower() == 'nob':
        db_list = ['nxrA','nxrB']
    for g_id in g_ids:
        q_file = glob(f'./genome_protein_files/{g_id}*')[0]
        for db_name in db_list:
            db = join('./curated_genes', db_name+'.faa')
            ofile = join(odir,'%s_%s' % (g_id,db_name))
            cmd = f'blastp -query {q_file} -outfmt 6 -max_hsps 1 -evalue 1e-3 -db {db} > {ofile}'
            check_call(cmd,shell=True)
        
# maually assigned
a = {'GCA_000341545.2': {'K00371': ['WP_018047899.1'], 'K00370': ['WP_042251421.1']},
     # for this, we need to manually add WP_018047899.1 into the sequenceing project.
     'GCA_000297255.1': {'K00370': ['CCF84486.1'], 'K00371': ['CCF85658.1']},
     'GCA_001597285.1': {'K10944': ['AMR68691.1']},   # haoC:AMR66258.1
     # 'GCA_003031105.1':{'K00370':['003031105_01085'],'k00371':['003031105_01084']},
     'GCA_000967305.2': {}  # haoA: ALO39610.1 , amoA: ALO39035.1
     }

sname2ko2locus.update(a)
for key in list(failed_g.keys()):
    failed_g_ids = failed_g.pop(key)
    for failed_g_id in failed_g_ids:
        if not judge(sname2ko2locus[failed_g_id], g_df.loc[failed_g_id, 'type']):
            failed_g[key].append(failed_g_id)

# maually assigned
target_g_id = 'GCA_001597285.1'
for db in glob(join('curated_genes', '*.phr')):
    db = db.replace('.phr', '')
    q_file = glob(f'./genome_protein_files/{target_g_id}*')[0]
    cmd = f'blastp -query {q_file} -outfmt 6 -max_hsps 1 -evalue 1e-3 -db {db} '
