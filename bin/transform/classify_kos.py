"""
Transform KO id like K00928\nK13283 into a table with some columns could help use to classify them
"""
import string
from collections import defaultdict

from bioservices import KEGG

kegg = KEGG()

target_dir = "/home-user/thliao/data/protein_db/kegg"


def ko_classified_br(ko_set):
    br2kos = defaultdict(list)
    mapping_file = f"{target_dir}/ko/ko_brite.list"
    for row in open(mapping_file):
        row = row.strip('\n')
        rows = row.split('\t')
        cur_ko = rows[0].split(':')[-1]

        if cur_ko in ko_set:
            br2kos[rows[1].split(':')[-1]].append(cur_ko)
    br2kos.pop('ko00001')  # the top of Brite hierarchy, useless...
    return br2kos


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

    infos = {}
    for br in iter_br:
        if not br.startswith('br'):
            new_br = f"br:{br}"
        hier_infos = kegg.http_get(f"get/{new_br}", frmt="txt").split('\n')
        br_name = hier_infos[0].split('\t')[-1]
        if '+' in br_name:
            br_name = suppl_dict[br]
        if isinstance(kos, list):
            iter_ko = kos
        elif isinstance(kos, dict):
            iter_ko = kos[br]

        for ko in iter_ko:
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
                    infos[ko] = hier_dict
                    break
    return infos
