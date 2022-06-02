#################################################################################
#### get relative ko with diamond output against kegg database.
####
####
#################################################################################

from bioservices.kegg import KEGG
from collections import defaultdict
import pandas as pd
from tqdm import tqdm
import click
import os
from os.path import exists, join
import pickle


def parse_id(ID, max_try=10):
    info_str = 0
    count_ = 0
    return_dict = {}
    while count_ <= max_try:
        info_str = kegg.get(ID)
        # if failed, it will return 400 or 403. still an int
        if isinstance(info_str, str):
            break
        count_ += 1
    if not isinstance(info_str, str):
        if '+' not in ID:
            return
        for sub_id in ID.split('+'):
            get_dict = parse_id(sub_id)
            if get_dict is None:
                continue
            return_dict.update(get_dict)
        return return_dict
    info_dict_list = [kegg.parse('ENTRY ' + each_str)
                      for each_str in info_str.split('\nENTRY ')
                      if each_str]
    for info_dict in info_dict_list:
        if not isinstance(info_dict, dict):
            print(info_dict)
            continue
        source_organism = info_dict.get("ORGANISM", 'unknown')
        entry = info_dict.get('ENTRY', 'unknown').split(' ')[0]
        if entry.startswith('ENTRY'):
            entry = [_
                     for _ in info_dict.get('ENTRY', 'unknown').split(' ')
                     if _][1]
        _cache = [ori for ori, _ in zip(ID.split('+'),
                                        ID.lower().split('+'))
                  if entry.lower() in _]
        if len(_cache) == 0:
            print(entry, ID)
            continue
        locus = _cache[0]
        Orthology = info_dict.get("ORTHOLOGY", None)
        if Orthology is not None:
            KO_id = ";".join(sorted(Orthology.keys()))
        else:
            KO_id = None
        NCBI_refID = info_dict.get("DBLINKS", {}).get("NCBI-ProteinID", None)
        uniprot_refID = info_dict.get("DBLINKS", {}).get("UniProt", None)

        AA_seq = info_dict.get("AASEQ", '').replace(' ', '')
        if not AA_seq:
            print('No Amine acid sequence detected. weird... for locus:', locus)
        return_dict[locus] = dict(ID=locus,
                                  ko=KO_id,
                                  ncbi_id=NCBI_refID,
                                  uniprot_refID=uniprot_refID,
                                  source_organism=source_organism,
                                  AA_seq=AA_seq)
    return return_dict


def get_KO_info(ID, max_try=10):
    info_str = 0
    count_ = 0
    return_dict = {}
    while count_ <= max_try:
        info_str = kegg.get(ID)
        # if failed, it will return 400 or 403. still an int
        if isinstance(info_str, str):
            break
        count_ += 1
    if not isinstance(info_str, str):
        if '+' not in ID:
            return
        for sub_id in ID.split('+'):
            get_dict = get_KO_info(sub_id)
            return_dict.update(get_dict)
        return return_dict
    info_dict_list = [kegg.parse('ENTRY ' + each_str)
                      for each_str in info_str.split('\nENTRY ')
                      if each_str]
    # make the first one entry startwith ENTRY instead of original locus.
    for info_dict in info_dict_list:
        if not isinstance(info_dict, dict):
            print(info_dict)
            continue
        entry = info_dict.get('ENTRY', 'unknown').split(' ')[0]
        if entry.startswith('ENTRY'):
            entry = [_
                     for _ in info_dict.get('ENTRY', 'unknown').split(' ')
                     if _][1]
        _cache = [ori for ori, _ in zip(ID.split('+'),
                                        ID.lower().split('+'))
                  if entry.lower() in _]
        if len(_cache) == 0:
            print(entry, ID)
            continue
        ko = _cache[0]
        gene_name = ';'.join(info_dict.get('NAME', ['']))
        definition = info_dict.get('DEFINITION', '')
        reference_t = ''
        if "REFERENCE" in info_dict:
            reference_t = ';'.join([str(_dict.get('TITLE', ''))
                                    for _dict in info_dict.get('REFERENCE', {})])

        return_dict[ko] = dict(gene_name=gene_name,
                               definition=definition,
                               reference_t=reference_t)
    return return_dict


def pack_it_up(ko2info, locus2ko, locus2info):
    tqdm.write("pack the result into a big dataframe")
    df_list = []
    for locus, ko_list in tqdm(locus2ko.items()):
        for ko in set(ko_list):
            ko_info = ko2info.get(ko, None)
            if ko_info is None:
                continue
            locus_info = locus2info[locus]
            _sub2 = pd.DataFrame().from_dict({locus: ko_info}, orient='index')
            _sub1 = pd.DataFrame().from_dict({locus: locus_info}, orient='index')
            _df = _sub1.join(_sub2, lsuffix=1)
            df_list.append(_df)
    tqdm.write("start concatenating......")

    if len(df_list) == 1:
        total_df = df_list[0]
    else:
        total_df = pd.concat(df_list, axis=0, sort=True)
    return total_df


def batch_iter(iter, batch_size):
    # generating batch according batch_size
    iter = list(iter)
    n_iter = []
    batch_d = 0
    for batch_u in range(0, len(iter), batch_size):
        if batch_u != 0:
            n_iter.append(iter[batch_d:batch_u])
        batch_d = batch_u
    n_iter.append(iter[batch_d: len(iter) + 1])
    return n_iter


@click.command(
    help="This script mainly for annotate diamond output against kegg databse. For using this script, please use python3.5+ and first install the `requirements`.\n\n just simply use python3 thisscript.py -i input_tab -o output_name.tsv ")
@click.option("-i", "input_tab")
@click.option("-o", "output_tab")
@click.option("-test", "test", is_flag=True, default=False)
def main(input_tab, output_tab, test):
    tmp_dir = './tmp'
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(output_tab)),
                exist_ok=True)
    df = pd.read_csv(input_tab, sep='\t', header=None)

    if test:
        # if use test option. just subtract top 50 to process for saving time.
        random50 = pd.np.random.choice(df.index, 50)
        df = df.loc[random50, :]
    tqdm.write("Get all relative information of the subject locus... ...")
    unique_DBlocus = set(df.loc[:, 1].unique())
    pack10_up = batch_iter(unique_DBlocus, 10)
    null_ID = []
    # 10 times faster
    if not exists(join(tmp_dir, 'dblocus2info')):
        DBlocus2info = {}
        for joined_DBlocus in tqdm(pack10_up, ):
            # todo: use asyncio to improve the speed
            DBlocus_info = parse_id('+'.join(joined_DBlocus))
            if DBlocus_info is None:
                null_ID += joined_DBlocus
                continue
            DBlocus2info.update(DBlocus_info)
        pickle.dump(DBlocus2info, open(join(tmp_dir, 'dblocus2info'), 'wb'))
        pickle.dump(null_ID, open(join(tmp_dir, 'null_ID'), 'wb'))
    else:
        DBlocus2info = pickle.load(open(join(tmp_dir, 'dblocus2info'), 'rb'))
        null_ID = pickle.load(open(join(tmp_dir, 'null_ID'), 'rb'))

    locus2info = {}
    locus2ko = defaultdict(list)
    ko2locus = defaultdict(list)
    for rid, row in df.iterrows():
        locus = row[0]
        DBlocus = row[1]
        if locus in locus2info:
            continue
        record = DBlocus2info.get(DBlocus, None)
        if record is None:
            continue
        if record.get("ko") is None:
            continue

        locus2info[row[0]] = record
        ko_list = record['ko']
        ko_list = ko_list.split(';')
        # some locus may not assigned with ko
        for ko in ko_list:
            locus2ko[locus].append(ko)
            ko2locus[ko].append(locus)

    ########################################################
    tqdm.write("collect all KO id, start iterate all KO info")
    if not exists(join(tmp_dir, 'ko2info')):
        ko2info = {}
        ko_list = list(ko2locus.keys())
        pack10_up = batch_iter(ko_list, 10)
        for ko_list in tqdm(pack10_up):
            ko_info = get_KO_info('+'.join(ko_list))
            if ko_info is None:
                continue
            ko2info.update(ko_info)
        pickle.dump(ko2info, open(join(tmp_dir, 'ko2info'), 'wb'))
    else:
        ko2info = pickle.load(open(join(tmp_dir, 'ko2info'), 'rb'))
    locus_df = pack_it_up(ko2info, locus2ko, locus2info)
    locus_df = locus_df.reindex(columns=['locus_tag',
                                         'ko',
                                         'definition',
                                         'gene_name',
                                         'ncbi_id',
                                         'uniprot_refID',
                                         'source_organism',
                                         'ID',
                                         'AA_seq',
                                         'reference_t'])
    locus_df.to_csv(output_tab, sep='\t', index=1, index_label='locus_tag')
    with open(output_tab + '.null_ID', 'w') as f1:
        f1.write('\n'.join(null_ID))
    return locus_df


if __name__ == '__main__':
    kegg = KEGG()
    main()
