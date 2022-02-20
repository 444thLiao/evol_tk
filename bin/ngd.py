#! /env/python
"""
For implement some missing feature(but i need) for ncbi-genome-download
"""

import os
import time
from collections import defaultdict
from glob import glob
from os.path import *

import click
import ncbi_genome_download as ngd
from ete3 import NCBITaxa
from ncbi_genome_download import NgdConfig
from ncbi_genome_download.core import select_candidates
from tqdm import tqdm
import pandas as pd
import io

HOME = os.getenv("HOME")
db_dir = f"{HOME}/data/NCBI/"

metadata_files_dir = f"{HOME}/.cache/ncbi-genome-download/"


# genbank_bacteria_assembly_summary.txt
def from_name2ids(phylum_name,
                  dataset='genbank',
                  return_d2ids=False):
    """
    retrieve ids and metadata from genbank file
    :param phylum_name:
    :return:
    """
    phylum_names = [_ for _ in phylum_name.split(';') if _]
    # phylum_name = "Nitrospirae;"
    # phylum_tid = "40117"
    ncbi = NCBITaxa()

    p2tid = ncbi.get_name_translator(phylum_names)

    for _ in phylum_names:
        if not p2tid.get(_):
            print(f" '{_}'' not found. please check the name")
    tids = [p2tid.get(_, [None])[0]
            for _ in phylum_names
            if p2tid.get(_)]
    tid2name = {p2tid.get(_, [None])[0]: _
                for _ in phylum_names
                if p2tid.get(_)}

    domain2dids = defaultdict(list)
    descend_ids = []
    tid2dids = {}
    for tid in tids:
        lineages = ncbi.get_lineage(tid)
        ranks = ncbi.get_rank(lineages)
        ranks = {v: k for k, v in ranks.items()}
        names = ncbi.get_taxid_translator(lineages)
        domain = names[ranks['superkingdom']]

        _descend_ids = ncbi.get_descendant_taxa(tid,
                                                intermediate_nodes=True)
        tid2dids[tid2name[tid]] = len(_descend_ids)
        descend_ids += _descend_ids
        domain2dids[domain].extend(_descend_ids)
    print(f"in total, {len(descend_ids)} taxids were found. ")
    if return_d2ids:
        return domain2dids

    domain2aids = defaultdict(list)
    collect_info = []
    descend_ids = set(descend_ids)
    for domain, ids in domain2dids.items():
        d = domain.lower()
        metadata = join(metadata_files_dir,
                        f"{dataset}_{d}_assembly_summary.txt")
        tqdm.write(
            f'read {metadata} which last modified at : {time.ctime(os.path.getmtime(metadata))}')
        for row in tqdm(open(metadata)):
            if row.startswith("GC"):
                rows = row.split('\t')
                if int(rows[5]) in descend_ids:
                    collect_info.append(row)
                    domain2aids[d].append(rows[0])
    return domain2aids, collect_info


def id2domain_to_ids(ids_list):
    domain2aids = defaultdict(list)
    collect_info = []
    ids_list = set(ids_list)
    for d in ["bacteria", 'archaea']:
        metadata = join(metadata_files_dir,
                        f"genbank_{d}_assembly_summary.txt")
        tqdm.write(f'read {metadata}')
        for row in tqdm(open(metadata)):
            if row.startswith("GC"):
                rows = row.split('\t')
                if str(rows[0]) in ids_list:
                    collect_info.append(row)
                    domain2aids[d].append(rows[0])
    missing_ids = ids_list.difference(
        set([_.split('\t')[0] for _ in collect_info]))
    if missing_ids:
        tqdm.write(f'{len(missing_ids)} are missing in summary file.')
    return domain2aids, collect_info


def from_tid2ids(taxons, dataset='genbank'):
    ncbi = NCBITaxa()

    def desc_taxa(taxid):
        descendent_taxa = ncbi.get_descendant_taxa(
            taxid, intermediate_nodes=True)
        descendent_taxas = []

        for _taxid in descendent_taxa:
            descendent_taxa_names = ncbi.translate_to_names([_taxid])
            descendent_taxas.append((str(_taxid), descendent_taxa_names[0]))

        return descendent_taxas

    def get_aid_from_tid(all_taxas):
        tids = [_[0] for _ in all_taxas]
        aids = []
        for d in ["bacteria", 'archaea']:
            metadata = join(metadata_files_dir,
                            f"{dataset}_{d}_assembly_summary.txt")
            a = pd.read_csv(io.StringIO(
                '\n'.join(open(metadata).read().split('\n')[1:])), sep='\t')
            aids += list(a.loc[a['taxid'].astype(str).isin(tids),
                               '# assembly_accession'])
        return aids

    all_taxas = []
    for tid in taxons.split(';'):
        descendent_taxas = desc_taxa(int(tid))
        all_taxas += descendent_taxas

    aids = get_aid_from_tid(all_taxas)
    return aids

# cids,cinfo = from_name2ids("Verrucomicrobia")


def batch_iter(iter, batch_size):
    # generating batch according batch_size
    iter = list(iter)
    n_iter = []
    batch_d = 0
    for batch_u in range(0, len(iter), int(batch_size)):
        if batch_u != 0:
            n_iter.append(iter[batch_d:batch_u])
        batch_d = batch_u
    n_iter.append(iter[batch_d: len(iter) + 1])
    return n_iter


def check_not_down(formats, acc_set, domain, odir):
    old_d = acc_set[::]
    refseq_acc = [_ for _ in old_d if _.startswith('GCF')]
    genbank_acc = [_ for _ in old_d if _.startswith('GCA')]

    format2suffix = {'fasta': '*.fna.gz',
                     "protein-fasta": "*.faa.gz", }
    sub_aids_all = []
    for f in formats:
        if f in format2suffix:
            suffix = format2suffix[f]
            # check whether other kinds of files have been downloaded
            sub_aids = []
            for acc in refseq_acc:
                if not glob(join(odir, 'refseq', domain, acc, suffix)):
                    sub_aids.append(acc)
            for acc in genbank_acc:
                if not glob(join(odir, 'genbank', domain, acc, suffix)):
                    sub_aids.append(acc)
            sub_aids_all.extend(sub_aids)
        else:
            sub_aids_all.extend(refseq_acc+genbank_acc)
    return list(set(sub_aids_all))


def main(name=None,
         odir=None,
         taxons=None,
         formats='fasta',
         ids_list=None,
         size_of_batch=30,
         parallel=10,
         enable_check=True,
         section='genbank',
         group='bacteria',
         dry_run=False):
    formats = formats.split(',')
    if odir is None:
        odir = db_dir
    else:
        odir = realpath(odir)
    if enable_check:
        if ids_list:
            # should be assembly ID list
            domain2aids, cinfos = id2domain_to_ids(ids_list)
        elif name is not None:
            domain2aids, cinfos = from_name2ids(name, dataset=section)
        elif taxons is not None:
            domain2aids, cinfos = from_tid2ids(taxons)

        # filter with existing files
        downloaded_aids = []
        new_domain2aids = {}
        for d, aids in domain2aids.items():
            sub_aids = check_not_down(formats, aids, d, odir)
            new_domain2aids[d] = sub_aids
            downloaded_aids.extend(new_domain2aids[d])
            tqdm.write(
                f"domain: {d}, original number of ids: {len(aids)}, now ids: {len(new_domain2aids[d])} ")
    elif not enable_check and ids_list:
        # disable the check and give a list of ids_list
        downloaded_aids = ids_list[::]
    if dry_run:
        with open(f'{odir}/downloaded_aids.list', 'w') as f1:
            f1.write('\n'.join(downloaded_aids))
    _d = {
        "dry_run": dry_run,
        "section": section,
        "groups": group,
        "parallel": parallel,
        "output": odir,
        "file_formats": formats}
    tqdm.write(f'params is {_d}')
    for batch_aids in tqdm(batch_iter(downloaded_aids, size_of_batch)):
        ngd.download(**{"assembly_accessions": ','.join(batch_aids),
                        "dry_run": dry_run,
                        "use_cache":True, # to avoid it automatic download/update the summary file 
                        "section": section,
                        "parallel": parallel,
                        "output": odir,
                        "groups": group,  # if not assign this, it will take long time to iterate all groups
                        "file_formats": formats})
    with open(join(odir, 'metadata.csv'), 'w') as f1:
        f1.write('\n'.join(cinfos))


@click.command(help="It would split the list of assembly ids into batches. The size parameter would produce")
@click.option("-n", "name", help="input the phylum name or other. use ; to separate multiple ")
@click.option("-t", "taxons", help="input the taxon id. It will retrieve all the genomes desceding to the provided taxon; to separate multiple ")
@click.option("-F", "formats", help='Which formats to download (default: %(default)s).'
                                    'A comma-separated list of formats is also possible. For example: "fasta,assembly-report". '
                                    'Choose from: {choices}'.format(
                                        choices=NgdConfig.get_choices('file_formats')),
              default='fasta')
@click.option("-o", "odir", help=f"Create output hierarchy in specified folder (default: {NgdConfig.get_default('output')})",
              default=NgdConfig.get_default('output'))
@click.option("-size", "size_of_batch", help=f"The size of each batch.",
              default=20)
@click.option("-p", "parallel", help=f"Run N downloads in parallel (default: 10)",
              default=5)
@click.option("-id", "id_list", help=f" ',' separated assembly ids or a single file  ",
              default='')
@click.option("-C", "enable_check", help=f"use summary file or use the id input directly",
              is_flag=True, default=True)
@click.option("-s", "section", help=f"refseq or genbank",
              default='genbank')
@click.option("-g", "group", help=f"The NCBI taxonomic groups to download (default: bacteria).",
              default='bacteria')
@click.option("-dry_run", "dry_run", help=f"If given, only output the id needed to be downloaded",
              is_flag=True, default=False)
def cli(name, odir, taxons, formats, size_of_batch, parallel, enable_check, id_list, section, group, dry_run):
    if exists(id_list):
        ids_list = [_ for _ in open(id_list).read().split('\n') if _]
    elif type(id_list) == str and id_list:
        ids_list = id_list.split(',')
    else:
        ids_list = None

    main(name=name,
         odir=odir,
         taxons=taxons,
         formats=formats,
         ids_list=ids_list,
         size_of_batch=size_of_batch,
         parallel=parallel,
         enable_check=enable_check,
         section=section, group=group, dry_run=dry_run
         )


if __name__ == '__main__':
    cli()
