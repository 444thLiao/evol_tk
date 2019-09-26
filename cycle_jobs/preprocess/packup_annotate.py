
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

name2ko = {'amoA': 'K10944',
           'amoB': 'K10945',
           'amoC': 'K10946',
           'hao': 'K10535',
           'nxrA': 'K00370',
           'nxrB': 'K00371'}
ko2name = dict([(v, k) for k, v in name2ko.items()])
# validatation
all_g_ids = [_ for _ in open(ref_id, 'r').read().split('\n') if _]
g_df = pd.read_excel(genome_info, index_col=0)
assert not set(all_g_ids).difference(g_df.index)


sname2ko2locus = dict()
kid_list = [_ for _ in open(target_file, 'r').read().split('\n') if _]
for ofile in tqdm(glob(join(odir, '*.out'))):
    sname = basename(ofile).replace('.out', '')
    if sname in all_g_ids:
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
        # return all([bool(ko2locus.get(k)) for k in ['K10946','K10945','K10944','K10535']])
        return all([bool(ko2locus.get(k)) for k in ['K10944', 'K10945', 'K10946', 'K10535']])
    elif type_g == 'AOA':
        return all([bool(ko2locus.get(k)) for k in ['K10944', 'K10945', 'K10946']])


def reupdate_dict(failed_g, sname2ko2locus):
    # reupdate the dict which stodge failed_g
    failed_g = failed_g.copy()
    for key in list(failed_g.keys()):
        failed_g_ids = failed_g.pop(key)
        for failed_g_id in failed_g_ids:
            if not judge(sname2ko2locus[failed_g_id], g_df.loc[failed_g_id, 'type']):
                failed_g[key].append(failed_g_id)
    return failed_g


og_tsv = './genome_protein_files/OrthoFinder/Results_Sep25/Orthogroups/Orthogroups.tsv'
og_df = pd.read_csv(og_tsv, sep='\t', low_memory=False, index_col=0)


def use_og_reannoate_(failed_g, sname2ko2locus):
    failed_g = failed_g.copy()
    sname2ko2locus = sname2ko2locus.copy()
    for key in list(failed_g.keys()):
        failed_g_ids = failed_g.pop(key)
        for failed_g_id in failed_g_ids:
            ko2locus = sname2ko2locus[failed_g_id]
            for ko, locus_l in ko2locus.items():
                if not locus_l:
                    backup_id = og_df.loc[ko2og[ko], failed_g_id]
                    if len(list(backup_id)) != 1:
                        continue
                    backup_id = list(backup_id)[0]
                    if pd.isna(backup_id):
                        continue
                    sname2ko2locus[failed_g_id][ko] = [backup_id]
            if not judge(sname2ko2locus[failed_g_id], g_df.loc[failed_g_id, 'type']):
                failed_g[key].append(failed_g_id)
    return failed_g, sname2ko2locus


def update_ko2og(sname2ko2locus, failed_g=[]):
    # auto implement first?
    # use Orthogroups to reimplement it
    if isinstance(failed_g, dict):
        id_list = [v for vlist in failed_g.values() for v in vlist]
    else:
        id_list = []
    sub_df = og_df.reindex(columns=all_g_ids)
    ko2og = defaultdict(list)
    ko2og2names = defaultdict(lambda: defaultdict(list))
    for g_id in sub_df.columns:
        if g_id in id_list:
            continue
        ko2locus = sname2ko2locus[g_id]
        for ko, locus_l in ko2locus.items():

            if locus_l:
                for l in locus_l:
                    ogs = sub_df.index[sub_df.loc[:, g_id].str.contains(
                        l, regex=False).fillna(False)]
                    if not list(ogs):
                        print(g_id, ko, locus_l)
                        break
                    og = list(ogs)[0]
                    if og not in ko2og[ko]:
                        ko2og[ko].append(og)
                    ko2og2names[ko][og].append((g_id, l))
    return ko2og, ko2og2names


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

ko2og, ko2og2names = update_ko2og(sname2ko2locus, failed_g)
failed_g, sname2ko2locus = use_og_reannoate_(failed_g, sname2ko2locus)


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
        q_file = f'./genome_protein_files/{g_id}.faa'
        ko2locus = sname2ko2locus[g_id]
        missing_name = [ko2name.get(k, '')
                        for k, l_list in ko2locus.items() if not l_list]
        for db_name in db_list:
            if db_name not in missing_name:
                continue
            db = join('./curated_genes', db_name+'.faa')
            ofile = join(odir, '%s_%s' % (g_id, db_name))
            cmd = f'blastp -query {q_file} -outfmt 6 -max_hsps 1 -evalue 1e-3 -db {db} > {ofile}'
            if not exists(ofile):
                check_call(cmd, shell=True)
            if os.path.getsize(ofile) != 0:
                _t = pd.read_csv(ofile, sep='\t', index_col=0, header=None)
                _t = _t.sort_values(10)
                if len(set(_t.index[:5])) == 1:
                    locus_ID = list(set(_t.index[:5]))[0]
                    print(type_g,
                          g_id,
                          db_name,
                          locus_ID)
                    sname2ko2locus[g_id][name2ko[db_name]] = [locus_ID]
                    manually_blast_r[g_id][name2ko[db_name]] = [locus_ID]
failed_g = reupdate_dict(failed_g, sname2ko2locus)
failed_g, sname2ko2locus = use_og_reannoate_(failed_g, sname2ko2locus)
# robust support (from paper/kofam scan)
# s = {'LS423452.1':{'K00371': ['SPS06750.1','SPS06997.1'], 'K00370': ['SPS06751.1','SPS06998.1']}}

# maually assigned
# a = {'GCA_000341545.2': {'K00371': ['WP_018047899.1'], 'K00370': ['WP_042251421.1']},
#      # for this, we need to manually add WP_018047899.1 into the sequenceing project.
#      'GCA_000297255.1': {'K00370': ['CCF84486.1'], 'K00371': ['CCF85658.1']},
#      'GCA_001597285.1': {'K10944': ['AMR68691.1'],},   # haoC:AMR66258.1
#      # 'GCA_003031105.1':{'K00370':['003031105_01085'],'k00371':['003031105_01084']},
#      'GCA_000967305.2': {'K00370':['ALO39018.1]}  # , amoA: ALO39035.1
#      }
a = {'GCA_002869885.2': {'K10535': ['THJ10221.1']},
     'GCA_000967305.2': {'K00370': ['ALO39018.1'], 'K10944': ['ALO39035.1']},
     'GCA_001597285.1': {'K100944': ['AMR68691.1']}
     }


sname2ko2locus.update(a)
failed_g = reupdate_dict(failed_g, sname2ko2locus)
failed_g, sname2ko2locus = use_og_reannoate_(failed_g, sname2ko2locus)
ko2og, ko2og2names = update_ko2og(sname2ko2locus, failed_g)

odir = 'json_dump'
os.makedirs(odir,exist_ok=True)
import json
with open(join(odir,'ko2og.json'),'w') as f1:
    json.dump(ko2og,f1)
with open(join(odir,'manually_blast_r.json'),'w') as f1:
    json.dump(manually_blast_r,f1)
with open(join(odir,'ko2og2names.json'),'w') as f1:
    json.dump(ko2og2names,f1)