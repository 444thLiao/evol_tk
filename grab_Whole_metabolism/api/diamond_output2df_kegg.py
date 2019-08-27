from bioservices.kegg import KEGG
from collections import defaultdict
import pandas as pd
from tqdm import tqdm
import click
import os


def parse_id(ID):
    info_dict = kegg.parse(kegg.get(ID))
    Orthology = info_dict.get("ORTHOLOGY", None)
    if Orthology is not None:
        KO_id = ";".join(sorted(Orthology.keys()))
    else:
        KO_id = None
    NCBI_refID = info_dict.get("DBLINKS", {}).get("NCBI-ProteinID", None)
    uniprot_refID = info_dict.get("DBLINKS", {}).get("UniProt", None)
    source_organism = info_dict["ORGANISM"]
    AA_seq = info_dict["AASEQ"].replace(' ', '')

    return_dict = dict(ko=KO_id,
                       ncbi_id=NCBI_refID,
                       uniprot_refID=uniprot_refID,
                       source_organism=source_organism,
                       AA_seq=AA_seq)
    return return_dict


def get_KO_info(ID):
    info_dict = kegg.parse(kegg.get(ID))
    gene_name = ';'.join(info_dict.get('NAME', ['']))
    definition = info_dict['DEFINITION']
    reference_t = ''
    if "REFERENCE" in info_dict:
        reference_t = ';'.join([_dict.get('title') for _dict in info_dict.get('REFERENCE')])

    return_dict = dict(gene_name=gene_name,
                       definition=definition,
                       reference_t=reference_t)
    return return_dict


def pack_it_up():
    pass


@click.command()
@click.option("-i", "input_tab")
def main(input_tab):
    df = pd.read_csv(input_tab, sep='\t', header=None)
    locus2dict = defaultdict(list)
    tqdm.write("Start to parse each ID into KO id set... ...")
    for rid, row in tqdm(df.iterrows(),
                         total=df.shape[0]):
        locus = row[0]
        ID = row[1]
        locus2dict[locus].append(parse_id(ID))

    ko2locus = defaultdict(list)
    locus2ko = defaultdict(list)
    for locus, info_dict_list in locus2dict.items():
        for info_dict in info_dict_list:
            ko_list = info_dict["KO_id"]
            if ko_list is not None:
                ko_list = ko_list.split(';')
            for ko in ko_list:
                ko2locus[ko].append(locus)
                locus2ko[locus].append(ko)
    # filter out ko??

    ########################################################
    tqdm.write("collect all KO id, start iterate all KO info")
    for ko, in ko2locus.keys():
        # todo: filter ko. not dealing with all ko
        ko_info = get_KO_info(ko)


if __name__ == '__main__':
    kegg = KEGG()
    main()
