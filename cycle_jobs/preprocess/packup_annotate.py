import json
from subprocess import check_call
from glob import glob
from os.path import join, dirname, exists, basename
import os
from tqdm import tqdm
import pandas as pd
from collections import defaultdict


target_file = './nitrification_list'
genome_info = './genome_info_full.xlsx'
odir = './annotated'

name2ko = {'amoA': 'K10944',
           'amoB': 'K10945',
           'amoC': 'K10946',
           'hao': 'K10535',
           'nxrA': 'K00370',
           'nxrB': 'K00371'}
ko2name = dict([(v, k) for k, v in name2ko.items()])

# validatation
g_df = pd.read_excel(genome_info, index_col=0)
g_df = g_df.loc[g_df.loc[:,'used']!='no',:]
all_g_ids = list(g_df.index)
# assert not set(all_g_ids).difference(g_df.index)


sname2ko2locus = dict()
kid_list = [_ for _ in open(target_file, 'r').read().split('\n') if _]
for ofile in tqdm(glob(join(odir, '*.out'))):
    sname = basename(ofile).replace('.out', '')
    if sname in all_g_ids:
        # only capture the reference genome annotation
        locus2ko = dict([(_.split('\t')[0], _.split('\t')[1])
                         for _ in open(ofile).read().split('\n')
                         if _])
        ko2locus = dict()
        for kid in kid_list:
            locus = [l for l, ks in locus2ko.items() if kid == ks]
            ko2locus[kid] = locus
        sname2ko2locus[sname] = ko2locus


# check
def judge(ko2locus, type_g):
    if type_g == 'comammox':
        return all([bool(v) for v in ko2locus.values()])
    elif type_g == 'NOB':
        # nxrB K00371 nxrA K00370
        return all([bool(ko2locus.get(k)) for k in ['K00371', 'K00370']])
    elif type_g == 'AOB':
        # amoA: K10944
        # amoB: K10945
        # amoC: K10946
        # hao: K10535
        return all([bool(ko2locus.get(k)) for k in ['K10944', 'K10945', 'K10946', 'K10535']])
    elif type_g == 'AOA':
        return all([bool(ko2locus.get(k)) for k in ['K10944', 'K10945', 'K10946']])

# iterative update part
def reupdate_dict(failed_g, sname2ko2locus):
    # reupdate the dict which stodge failed_g 
    # just update failed_g with sname2ko2locus
    failed_g = failed_g.copy()
    for key in list(failed_g.keys()):
        failed_g_ids = failed_g.pop(key)
        for failed_g_id in failed_g_ids:
            if not judge(sname2ko2locus[failed_g_id], g_df.loc[failed_g_id, 'type']):
                failed_g[key].append(failed_g_id)
    return failed_g

def use_og_reannoate_(failed_g, sname2ko2locus):
    # update both failed_g, sname2ko2locus with og table
    failed_g = failed_g.copy()
    sname2ko2locus = sname2ko2locus.copy()
    for key in list(failed_g.keys()):
        failed_g_ids = failed_g.pop(key)
        for failed_g_id in failed_g_ids:
            ko2locus = sname2ko2locus[failed_g_id]
            for ko, locus_l in ko2locus.items():
                if not locus_l:
                    backup_ids = og_df.loc[ko2og[ko], failed_g_id]
                    backup_ids = [_ for _ in list(backup_ids) if not pd.isna(_)]
                    if len(backup_ids) != 1:
                        # multiple choice or no choice will pass it
                        continue
                    backup_id = backup_ids[0]
                    sname2ko2locus[failed_g_id][ko] = [backup_id]
            if not judge(sname2ko2locus[failed_g_id], g_df.loc[failed_g_id, 'type']):
                failed_g[key].append(failed_g_id)
    
    return failed_g, sname2ko2locus

def update_ko2og(sname2ko2locus, failed_g=[]):
    # get ko2og with defined failed_g
    # would not use og info from the ids of failed_g
    if isinstance(failed_g, dict):
        id_list = [v for vlist in failed_g.values() for v in vlist]
    else:
        id_list = []
    sub_df = og_df.reindex(columns=all_g_ids)
    # only consider these reference genomes
    ko2og = defaultdict(list)
    ko2og2names = defaultdict(lambda: defaultdict(list))
    for g_id in sub_df.columns:
        if g_id in id_list:
            continue
        ko2locus = sname2ko2locus[g_id]
        for ko, locus_l in ko2locus.items():
            if locus_l:
                for l in locus_l:
                    ogs = sub_df.index[sub_df.loc[:, g_id].fillna(
                        '').str.contains(l, regex=False)]
                    ogs = list(ogs)
                    if not ogs:
                        print(g_id, ko, locus_l)
                        break
                    og = ogs[0]
                    if og not in ko2og[ko]:
                        ko2og[ko].append(og)
                    ko2og2names[ko][og].append((g_id, l))
    return ko2og, ko2og2names

og_tsv = './genome_protein_files/OrthoFinder/Results_Sep27/Orthogroups/Orthogroups.tsv'
og_tsv = './genome_protein_files_more/OrthoFinder/Results_Oct01/Orthogroups/Orthogroups.tsv'
og_df = pd.read_csv(og_tsv, sep='\t', low_memory=False, index_col=0)

def rename(x):
    if pd.isna(x):
        return x
    if ', ' not in x:
        return str(x).split(' ')[0]
    else:
        return ', '.join([str(_).split(' ')[0] for _ in x.split(', ')])
og_df = og_df.applymap(rename)

# first time to get ko2og, no failed_g, all info are comed from kofamscan
ko2og, ko2og2names = update_ko2og(sname2ko2locus)

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

# reannotate the sname2ko2locus with og table
# re get ko2og, because the failed_g may be implemented by og table.
failed_g, sname2ko2locus = use_og_reannoate_(failed_g, sname2ko2locus)
ko2og, ko2og2names = update_ko2og(sname2ko2locus)


# manuall blast
manually_blast_r = defaultdict(dict)
odir = './reannotate'
os.makedirs(odir, exist_ok=True)
for type_g, g_ids in failed_g.items():
    if type_g == 'comammox':
        db_list = ['amoA', 'amoB', 'amoC', 'hao',
                   'haoA', 'HAO'] + ['nxrA', 'nxrB']
    elif type_g == 'AOB':
        db_list = ['amoA', 'amoB', 'amoC', 'hao', 'haoA', 'HAO']
    elif type_g == 'AOA':
        db_list = ['amoA', 'amoB', 'amoC']
    elif type_g == 'NOB':
        db_list = ['nxrA', 'nxrB']
    for g_id in g_ids:
        if g_id == 'GCA_002869895.2':
            continue # manually removed this one
        q_file = f'./genome_protein_files/{g_id}.faa'
        ko2locus = sname2ko2locus[g_id]
        missing_name = [ko2name.get(k, '')
                        for k, l_list in ko2locus.items() if not l_list]
        for db_name in db_list:
            if db_name not in missing_name:
                continue
            db = join('./curated_genes', db_name+'.faa')
            ofile = join(odir, '%s_%s' % (g_id, db_name))
            cmd = f'blastp -query {q_file} -outfmt 6 -max_hsps 1 -evalue 1e-4 -db {db} > {ofile}'
            if not exists(ofile):
                check_call(cmd, shell=True)
            if os.path.getsize(ofile) != 0:
                _t = pd.read_csv(ofile, sep='\t', index_col=0, header=None)
                _t = _t.sort_values(10)
                if len(set(_t.index[:5])) == 1:
                    locus_ID = list(set(_t.index[:5]))[0]
                    print('successfully reannotated', type_g,
                          g_id,
                          db_name,
                          locus_ID)
                    sname2ko2locus[g_id][name2ko[db_name]] = [locus_ID]
                    manually_blast_r[g_id][name2ko[db_name]] = [locus_ID]

# after update the sname2ko2locus
# update the failed_g
failed_g = reupdate_dict(failed_g, sname2ko2locus)
# reannotated with og table
failed_g, sname2ko2locus = use_og_reannoate_(failed_g, sname2ko2locus)
# re get ko2og
ko2og, ko2og2names = update_ko2og(sname2ko2locus, failed_g)
# re defined failed_g
failed_g, sname2ko2locus = use_og_reannoate_(failed_g, sname2ko2locus)


# robust support (from paper/kofam scan)
# s = {'LS423452.1':{'K00371': ['NITFAB_2343','NITFAB_2595'], 'K00370': ['NITFAB_2344','NITFAB_2596']}}

# maually assigned
# GCA_000341545.2 for this, we need to manually add WP_018047899.1 into the sequenceing project.

a = {
     'GCA_000967305.2': {'K10944': ['ALO39035.1','ALO39034.1'],'K10535': ['ALO37589.1']},
     'GCA_001597285.1': {'K10944': ['AMR68691.1'],}, 
     'GCA_000007565.2': {'K10944':['AAN67037.1','AAN67038.1'],},
     #'GCA_002869925.2':{'K10945':['THI84063.1']}   # fake one....
     }
manually_blast_r.update(a)

sname2ko2locus.update(a)
failed_g = reupdate_dict(failed_g, sname2ko2locus)
failed_g, sname2ko2locus = use_og_reannoate_(failed_g, sname2ko2locus)
ko2og, ko2og2names = update_ko2og(sname2ko2locus, failed_g)

odir = 'json_dump_v2'
os.makedirs(odir, exist_ok=True)
with open(join(odir, 'ko2og.json'), 'w') as f1:
    json.dump(ko2og, f1)
with open(join(odir, 'manually_blast_r.json'), 'w') as f1:
    json.dump(manually_blast_r, f1)
with open(join(odir, 'ko2og2names.json'), 'w') as f1:
    json.dump(ko2og2names, f1)
