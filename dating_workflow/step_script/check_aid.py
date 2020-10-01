"""
It main return counts of the existence of each steps of given ID list.
Including following steps:
1. genomes
2. proteins (generated by prokka) / or the results of prokka
3. bac120 annotations
4. cog25 annotations

"""
import os
from glob import glob
from os.path import *

import click
import pandas as pd

HOME = os.getenv("HOME")
db_path = f"{HOME}/data/NCBI/genbank/"
proteins_path = f"{HOME}/data/NCBI/modified_data/genome_protein_files"

bac120_anno_path = f"{HOME}/data/NCBI/modified_data/bac120_annotate"
cog25_anno_path = f"{HOME}/data/NCBI/modified_data/cog25_annotate"


def _check(p):
    if not p:
        return False
    else:
        if getsize(p[0]) or isdir(p[0]):
            return True
        else:
            return False


def check_genome(aid):
    f = glob(join(db_path, "*", aid))
    return _check(f)


def check_protein(aid):
    f = glob(join(proteins_path, f"{aid}.faa"))
    return _check(f)


def check_prokka(aid):
    f = glob(join(proteins_path, '../prokka_o', f"{aid}"))
    return _check(f)


def check_bac120(aid):
    f1 = glob(join(bac120_anno_path, 'PFAM', f"{aid}.out"))
    f2 = glob(join(bac120_anno_path, 'TIGRFAM', f"{aid}.out"))
    if _check(f1) and _check(f2):
        return True
    else:
        return False


def check_cog25(aid):
    f = glob(join(cog25_anno_path, f"{aid}.out"))
    return _check(f)


def main(genome_list):
    # test the existence of genome
    all_ids = open(genome_list).read().split('\n')
    all_ids = [_ for _ in all_ids if _]
    df = pd.DataFrame(index=['counts'],
                      data=0,
                      columns=["genomes", 'protein', 'prokka',
                               'bac120', 'cog25',
                               # r55, 16S...
                               ])
    no_genomes_ids = []
    no_protein_ids = []
    no_prokka_ids = []
    no_bac120_ids = []
    no_cog25_ids = []
    for aid in all_ids:
        if not check_genome(aid):
            no_genomes_ids.append(aid)
        if not check_protein(aid):
            no_protein_ids.append(aid)
        if not check_prokka(aid):
            no_prokka_ids.append(aid)
        if not check_bac120(aid):
            no_bac120_ids.append(aid)
        if not check_cog25(aid):
            no_cog25_ids.append(aid)
    stats = [no_genomes_ids,
             no_prokka_ids,
             no_prokka_ids,
             no_bac120_ids,
             no_cog25_ids]
    df.loc['counts', :] = list(map(len, stats))
    return df, stats


@click.command()
@click.option("-gl", "genome_list", default=None,
              help="It will read 'selected_genomes.txt', please prepare the file, or indicate the alternative name or path. It could be None. If you provided, you could use it to subset the aln sequences by indicate names.")
@click.option("-o", "ofile", default='./check.tab')
def cli(genome_list, ofile):
    df, stats = main(genome_list)
    df.to_csv(ofile, sep='\t', index_label="assembly ID", index=1)


if __name__ == '__main__':
    cli()
