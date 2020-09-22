import io
import os
from collections import defaultdict
from os.path import *

import click
import pandas as pd
from ete3 import NCBITaxa

ncbi = NCBITaxa()

HOME = os.getenv("HOME")
base_file = f"{HOME}/.cache/ncbi-genome-download/genbank_bacteria_assembly_summary.txt"
a = pd.read_csv(io.StringIO('\n'.join(open(base_file).read().split('\n')[1:])), sep='\t')

#
# taxid = '914'
# tax_name = ncbi.get_taxid_translator([taxid])


def process_path(path):
    if not '/' in path:
        path = './' + path
    if path.startswith('~'):
        path = expanduser(path)
    if path.startswith('.'):
        path = abspath(path)
    return path


def desc_taxa(taxid):
    descendent_taxa = ncbi.get_descendant_taxa(taxid,intermediate_nodes=True)
    descendent_taxas = []

    for _taxid in descendent_taxa:
        descendent_taxa_names = ncbi.translate_to_names([_taxid])
        descendent_taxas.append((str(_taxid), descendent_taxa_names[0]))

    return descendent_taxas


def get_aid_from_tid(all_taxas):
    tids = [_[0] for _ in all_taxas]
    aids = a.loc[a['taxid'].astype(str).isin(tids), :]
    return aids


@click.command()
@click.option("-t", "taxons", default="914")
@click.option("-n", "names", default="Nitrosomonas")
@click.option("-o", 'odir', default='./')
def cli(taxons, names, odir):
    odir = process_path(odir)
    if not exists(odir):
        os.makedirs(odir)
    tids = taxons.split(',')
    names = names.split(',')
    name2tid = dict(zip(names, tids))

    all_taxas = []
    for n, tid in name2tid.items():
        descendent_taxas = desc_taxa(int(tid))
        all_taxas += descendent_taxas

    aids_sub_df = get_aid_from_tid(all_taxas)
    # stats
    tid2count = defaultdict(int)
    for tid in aids_sub_df['taxid']:
        lineages = ncbi.get_lineage(tid)
        lineages = [str(_) for _ in lineages]
        for _tid in name2tid.values():
            if _tid in lineages:
                tid2count[_tid] += 1

    aids = aids_sub_df.iloc[:, 0]
    with open(join(odir, 'request_aids.list'), 'w') as f1:
        f1.write('\n'.join(aids))
    with open(join(odir, 'request.log'), 'w') as f1:
        f1.write("request tid and names are listed below: " + '\n')
        f1.write("name\ttid\taids count\n")
        for n, t in name2tid.items():
            f1.write('\t'.join([n, t, str(tid2count[t])]) + '\n')
    # taxons = '206379,914,1293497,40117'
    # names = 'Nitosomonadaceae,Nitrosomonas,Nitrospinae,Nitrospirae'


if __name__ == '__main__':
    cli()
