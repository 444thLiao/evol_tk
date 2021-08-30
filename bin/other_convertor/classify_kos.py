"""
Transform KO id like K00928\nK13283 into a table with some columns could help use to classify them

from bin.other_convertor.classify_kos import get_md_infos

"""
import string
from collections import defaultdict

import pandas as pd
from bioservices import KEGG
from tqdm import tqdm
import os

kegg = KEGG()

HOME = os.getenv('HOME')
target_dir = f"{HOME}/data/protein_db/kegg"


def ko_classified_br(ko_set):
    br2kos = defaultdict(list)
    mapping_file = f"{target_dir}/ko/ko_brite.list"
    for row in open(mapping_file):
        row = row.strip('\n')
        rows = row.split('\t')
        cur_ko = rows[0].split(':')[-1]

        if cur_ko in ko_set:
            br2kos[rows[1].split(':')[-1]].append(cur_ko)
    # br2kos.pop('ko00001')  # the top of Brite hierarchy, useless...
    return br2kos


def ko_classified_module(ko_set):
    md2kos = defaultdict(list)
    mapping_file = f"{target_dir}/ko/ko_module.list"
    for row in open(mapping_file):
        row = row.strip('\n')
        rows = row.split('\t')
        cur_ko = rows[0].split(':')[-1]

        if cur_ko in ko_set:
            md2kos[rows[1].split(':')[-1]].append(cur_ko)
    # br2kos.pop('ko00001')  # the top of Brite hierarchy, useless...
    return md2kos


def get_md_infos(mds):
    if isinstance(mds, str):
        mds = set([mds])
    else:
        mds = set(mds)
    md2info = {}

    for md in tqdm(mds):
        text = 404
        while type(text) == int:
            text = kegg.get(md)
        text = text.split('\n')
        name = [_.split('  ')[-1] for _ in text if _.startswith('NAME')][-1]
        all_kos = [_.split('  ')[-1] for _ in text if _.startswith('DEFINITION')][-1]
        # all_kos = [ko for e in all_kos.strip().split(' ') for ko in e.split('+')]
        md2info[md] = {'name': name,
                       'all ko': all_kos}
    return md2info

def ko2pathway(kos,max_try=50):
    if type(kos) == str:
        kos = [kos]

    ko2p = {}

    for ko in tqdm(kos):
        pathway_collect = []
        c = False
        text = 404
        _count = 0
        while type(text) == int:
            text = kegg.get(ko)
            _count +=1
            if _count>=max_try:
                break
        if type(text)==int:
            continue
        for row in text.split('\n'):
            if row.startswith("PATHWAY"):
                c = True
            if row.startswith("MODULE") or row.startswith("BRITE") or row.startswith("DISEASE"):
                c = False
            if c:
                pathway_collect.append(row.split('PATHWAY')[-1].strip())
        pathway_collect = {_.split(' ')[0].strip():_.partition(' ')[-1].strip() for _ in pathway_collect}
        ko2p[ko] = pathway_collect
    return ko2p


def get_ko_infos(kos):
    if isinstance(kos, str):
        kos = set([kos])
    else:
        kos = set(kos)
    kegg_info_tab = f"{target_dir}/ko_info.tab"
    ko2info_dict = {row.split('\t')[0].split(':')[-1]: row.strip('\n').split('\t')[1]
                    for row in open(kegg_info_tab)
                    if row.split('\t')[0].split(':')[-1] in kos}
    return ko2info_dict


suppl_dict = {'ko02000': 'Transporters',
              "ko03029": "Mitochondrial Biogenesis"}


def get_br_info(br, kos=None):
    if isinstance(br, str):
        iter_br = [br]
    elif isinstance(br, dict):
        iter_br = list(br)
        kos = br.copy()
    else:
        iter_br = list(br)

    infos = []
    for br in iter_br:
        if not br.startswith('br'):
            new_br = f"br:{br}"
        else:
            continue
        hier_infos = str(kegg.http_get(f"get/{new_br}", frmt="txt")).split('\n')
        if not hier_infos:
            print(br)
        br_name = hier_infos[0].split('\t')[-1]
        if '+' in br_name:
            br_name = suppl_dict.get(br, 'MISSING')
        if isinstance(kos, list):
            iter_ko = kos
        elif isinstance(kos, dict):
            iter_ko = kos[br]

        for ko in tqdm(iter_ko):
            hier_dict = {'top': br_name, 'top_br': br}
            level_index = string.ascii_letters[26:]  # ABCDEF...
            for row in hier_infos[1:]:
                if not row:
                    continue
                if row[0] in level_index:
                    hier_dict[row[0]] = row[1:].strip()
                if ko in row:
                    for _ in level_index[level_index.index(row[0]):]:
                        if _ in hier_dict:
                            hier_dict.pop(_)
                    _cache = pd.DataFrame.from_dict({ko: hier_dict}, orient='index')
                    infos.append(_cache)
                    # infos[ko] = hier_dict
                    # break
                    # if break, only the first hierarchy would be retrieved
    return infos


