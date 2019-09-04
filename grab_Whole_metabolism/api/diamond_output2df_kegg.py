from bioservices.kegg import KEGG
from collections import defaultdict, Counter
import pandas as pd
from tqdm import tqdm
import click
import os
from os.path import exists, join
import pickle
import time


def parse_id(ID, max_try=500):
    info_str = 0
    count_ = 0
    while count_ <= max_try:
        info_str = kegg.get(ID)
        # if failed, it will return 400 or 403. still an int
        if isinstance(info_str, str):
            break
        count_ += 1
        time.sleep(0.01)
    if not isinstance(info_str, str):
        return
    assert info_str.count('ENTRY') == len(ID.split('+'))
    info_dict_list = [kegg.parse('ENTRY' + each_str)
                      for each_str in info_str.split('ENTRY')
                      if each_str]
    return_dict = {}
    for locus, info_dict in zip(ID.split('+'),
                                info_dict_list):
        source_organism = info_dict.get("ORGANISM", 'unknown')
        Orthology = info_dict.get("ORTHOLOGY", None)
        if Orthology is not None:
            KO_id = ";".join(sorted(Orthology.keys()))
        else:
            KO_id = None
        NCBI_refID = info_dict.get("DBLINKS", {}).get("NCBI-ProteinID", None)
        uniprot_refID = info_dict.get("DBLINKS", {}).get("UniProt", None)
        source_organism = info_dict.get("ORGANISM", 'unknown')
        AA_seq = info_dict.get("AASEQ", '').replace(' ', '')
        if not AA_seq:
            print('No Amine acid sequence detected. weird... for ID:', ID)
        return_dict[locus] = dict(ID=ID,
                                  ko=KO_id,
                                  ncbi_id=NCBI_refID,
                                  uniprot_refID=uniprot_refID,
                                  source_organism=source_organism,
                                  AA_seq=AA_seq)
    return return_dict


def get_KO_info(ID, max_try=500):
    info_str = 0
    count_ = 0
    while count_ <= max_try:
        info_str = kegg.get(ID)
        # if failed, it will return 400 or 403. still an int
        if isinstance(info_str, str):
            break
        count_ += 1
        time.sleep(0.01)
    if not isinstance(info_str, str):
        return
    assert info_str.count('ENTRY') == len(ID.split('+'))
    info_dict_list = [kegg.parse('ENTRY' + each_str)
                      for each_str in info_str.split('ENTRY')
                      if each_str]
    return_dict = {}
    for ko, info_dict in zip(ID.split('+'),
                             info_dict_list):
        gene_name = ';'.join(info_str.get('NAME', ['']))
        definition = info_str.get('DEFINITION', '')
        reference_t = ''
        if "REFERENCE" in info_str:
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
        for ko in ko_list:
            ko_info = ko2info.get(ko, None)
            if ko_info is None:
                continue
            locus_info_list = locus2info[locus]
            for locus_info in locus_info_list:
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
@click.option("-no_highest", "get_highest", help='default is getting highest one for annotating.', is_flag=True, default=True)
@click.option("-drop_dup_ko", "drop_dup_ko", help='drop duplicated ko for single query subject.', is_flag=True, default=False)
@click.option("-test", "test", is_flag=True, default=False)
def main(input_tab, output_tab, get_highest, drop_dup_ko, test):
    tmp_dir = './tmp'
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(output_tab)),
                exist_ok=True)
    df = pd.read_csv(input_tab, sep='\t', header=None)
    if get_highest:
        df = df.sort_values([0, 10], ascending=True)
        df = df.drop_duplicates(0)
        # get smallest e-value one, and drop others
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

    locus2info = defaultdict(list)
    for rid, row in df.iterrows():
        locus2info[row[0]].append(DBlocus2info[row[1]])

    locus2ko = defaultdict(list)
    ko2locus = defaultdict(list)
    for locus, info_dict_list in locus2info.items():
        for info_dict in info_dict_list:
            ko_list = info_dict["ko"]
            if ko_list is not None:
                ko_list = ko_list.split(';')
            else:
                continue
            # some locus may not assigned with ko
            for ko in ko_list:
                locus2ko[locus].append(ko)
                ko2locus[ko].append(locus)

    if drop_dup_ko:
        # update locus2ko and ko2locus
        ko2locus = defaultdict(list)
        _locus2ko = dict()
        for locus, ko_list in locus2ko.items():
            # choose only one ko for each locus.
            num_ko = len(ko_list)
            freq_ko = {k: v / num_ko
                       for k, v in Counter(ko_list).items()}
            lg_60 = [k for k, v in freq_ko.items() if v >= 0.6]
            if len(lg_60) == 1:
                ko2locus[lg_60[0]].append(locus)
                _locus2ko[locus] = [lg_60[0]]
            elif len(lg_60) >= 2:
                # impossible....
                pass
            else:
                # no larger than 60%
                print("no large than 60, locus : {0}".format(locus))
        locus2ko = _locus2ko

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
    locus_df.to_csv(output_tab, sep='\t', index=1, index_label='locus_tag')
    with open(output_tab + '.null_ID', 'w') as f1:
        f1.write('\n'.join(null_ID))
    return locus_df


if __name__ == '__main__':
    kegg = KEGG()
    main()
