"""
It contain two kinds of functions.

1. convert the summary tab into a taxonomy tab
2. convert a assembly ID into taxid throught Entrez/api.

"""
import os

import pandas as pd
from tqdm import tqdm

from bin.ncbi_convertor import tax2tax_info
from api_tools.IO_for.read import read_summary
HOME = os.getenv("HOME")
db_dir = f"{HOME}/data/NCBI/"

default_infile = f"{HOME}/.cache/ncbi-genome-download/genbank_bacteria_assembly_summary.txt"

default_taxon_tab = f"{HOME}/.cache/ncbi-genome-download/taxonomy.tab"


def file2taxon_tab(infile):
    # infile = "/home-user/thliao/.cache/ncbi-genome-download/genbank_bacteria_assembly_summary.txt"
    metadata_df = read_summary(infile)
    tax_df = pd.read_csv(default_taxon_tab, sep='\t', index_col=0)
    missing_gids = [_ for _ in metadata_df.index if _.split('.')[0] not in tax_df.index]
    print(len(missing_gids))
    sub_df = metadata_df.loc[missing_gids]
    failed_tids = []
    aid2taxon_info = {}
    for aid, row in tqdm(sub_df.iterrows(), total=sub_df.shape[0]):
        try:
            taxon_dict = tax2tax_info(int(row['taxid']))
            aid2taxon_info[aid] = taxon_dict
        except:
            failed_tids.append(int(row['taxid']))
    print(f"failed taxids: {len(failed_tids)}")
    return aid2taxon_info

def rewrite_existing_tab(aid2taxon_info):
    _df = pd.DataFrame.from_dict(aid2taxon_info, orient='index')
    _df.index = [_.split('.')[0] for _ in _df.index]
    now_num = _df.shape[0]
    df2 = pd.read_csv(default_taxon_tab, sep='\t', index_col=0)
    df = df2.reindex(set(_df.index).union(df2.index))
    df.update(df2)
    df.update(_df)
    new_num = df.shape[0]
    print(f"latest tab contain {now_num} assembly IDs, now, merged tab contain {new_num} assembly IDs")
    df = df.sort_index()
    df.to_csv(default_taxon_tab,sep='\t', index=1, index_label=df.index.name)


# def aid2taxon(id_list, redo=False):
#     convertor = NCBI_convertor(id_list, "assembly")
#     suffix = 'protein2GI'
#     convertor.check_cache(suffix=suffix, redo=redo)
#     id2taxon = convertor.get_taxons_from_tid()
#     return id2taxon


def cli(infile=default_infile):
    aid2taxon_info = file2taxon_tab(infile)
    rewrite_existing_tab(aid2taxon_info=aid2taxon_info)


if __name__ == '__main__':
    import sys
    for f in sys.argv[1:]:
        cli(f)
