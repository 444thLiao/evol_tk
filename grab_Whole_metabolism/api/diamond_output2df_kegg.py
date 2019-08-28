from bioservices.kegg import KEGG
from collections import defaultdict, Counter
import pandas as pd
from tqdm import tqdm
import click
import os


def parse_id(ID, max_try=5):
    info_dict = 0
    count_ = 0
    while isinstance(info_dict, int):
        info_dict = kegg.parse(kegg.get(ID))
        # if failed, it will return 400 or 403. still an int
        count_ += 1
        if count_ >= max_try:
            return ID
    Orthology = info_dict.get("ORTHOLOGY", None)
    if Orthology is not None:
        KO_id = ";".join(sorted(Orthology.keys()))
    else:
        KO_id = None
    NCBI_refID = info_dict.get("DBLINKS", {}).get("NCBI-ProteinID", None)
    uniprot_refID = info_dict.get("DBLINKS", {}).get("UniProt", None)
    source_organism = info_dict["ORGANISM"]
    AA_seq = info_dict["AASEQ"].replace(' ', '')

    return_dict = dict(ID=ID,
                       ko=KO_id,
                       ncbi_id=NCBI_refID,
                       uniprot_refID=uniprot_refID,
                       source_organism=source_organism,
                       AA_seq=AA_seq)
    return return_dict


def get_KO_info(ID, max_try=5):
    info_dict = 0
    count_ = 0
    while isinstance(info_dict, int):
        info_dict = kegg.parse(kegg.get(ID))
        # if failed, it will return 400 or 403. still an int
        count_ += 1
        if count_ >= max_try:
            return ID
    gene_name = ';'.join(info_dict.get('NAME', ['']))
    definition = info_dict['DEFINITION']
    reference_t = ''
    if "REFERENCE" in info_dict:
        reference_t = ';'.join([_dict.get('TITLE', '')
                                for _dict in info_dict.get('REFERENCE')])

    return_dict = dict(gene_name=gene_name,
                       definition=definition,
                       reference_t=reference_t)
    return return_dict


def pack_it_up(ko2info, locus2ko, locus2info):
    total_df = pd.DataFrame()
    for locus, ko_list in locus2ko.items():
        for ko in ko_list:
            ko_info = ko2info[ko]
            locus_info_list = locus2info[locus]
            for locus_info in locus_info_list:
                _sub2 = pd.DataFrame().from_dict({locus: ko_info}, orient='index')
                _sub1 = pd.DataFrame().from_dict({locus: locus_info}, orient='index')
                _df = _sub1.join(_sub2, lsuffix=1)
                total_df = total_df.append(_df)
    return total_df


@click.command()
@click.option("-i", "input_tab")
@click.option("-o", "output_tab")
@click.option("-no_highest", "get_highest", is_flag=True, default=True)
@click.option("-drop_dup_ko", "drop_dup_ko", is_flag=True, default=False)
def main(input_tab, output_tab, get_highest, drop_dup_ko):
    os.makedirs(os.path.dirname(os.path.abspath(output_tab)),
                exist_ok=True)
    df = pd.read_csv(input_tab, sep='\t', header=None)
    if get_highest:
        df = df.sort_values([0, 10], ascending=True)
        df = df.drop_duplicates(0)
        # get smallest e-value one, and drop others
    df = df.iloc[:100, :]
    tqdm.write("Get all relative information of the subject locus... ...")
    unique_DBlocus = set(df.loc[:, 1].unique())
    DBlocus2info = {}
    null_ID = []
    for DBlocus in tqdm(unique_DBlocus,
                        total=len(unique_DBlocus)):
        # todo: use asyncio to improve the speed
        DBlocus_info = parse_id(DBlocus)
        if not isinstance(DBlocus_info, dict):
            null_ID.append(DBlocus_info)
        else:
            DBlocus2info[DBlocus] = DBlocus_info
    locus2info = {row[0]: DBlocus2info[row[1]]
                  for rid, row in df.iterrows()}

    locus2ko = defaultdict(list)
    ko2locus = defaultdict(list)
    for locus, info_dict_list in locus2info.items():
        for info_dict in info_dict_list:
            ko_list = info_dict["ko"]
            if ko_list is not None:
                ko_list = ko_list.split(';')
            for ko in ko_list:
                locus2ko[locus].append(ko)
                ko2locus[ko].append(locus)
    if drop_dup_ko:
        ko2locus = defaultdict(list)
        _locus2ko = dict()
        for locus, ko_list in locus2ko.items():
            # choose only one ko for each locus.
            num_ko = len(ko_list)
            freq_ko = {k: v / num_ko for k, v in Counter(ko_list).items()}
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
    ko2info = {}
    for ko, locus_list in tqdm(ko2locus.items(),
                               total=len(ko2locus)):
        ko_info = get_KO_info(ko)
        ko2info[ko] = ko_info
    locus_df = pack_it_up(ko2info, locus2ko, locus2info)
    locus_df.to_csv(output_tab, sep='\t', index=1, index_label='locus_tag')
    return locus_df


if __name__ == '__main__':
    kegg = KEGG()
    main()
