"""
For implement some missing feature(but i need) for ncbi-genome-download
"""

from os.path import *

import click
import ncbi_genome_download as ngd
from ete3 import NCBITaxa
from ncbi_genome_download import NgdConfig
from tqdm import tqdm
from collections import defaultdict
import os
HOME = os.getenv("HOME")
db_dir = f"{HOME}/data/NCBI/"

metadata_files_dir = f"{HOME}/.cache/ncbi-genome-download/"
# genbank_bacteria_assembly_summary.txt
def from_name2ids(phylum_name,
                  return_d2ids = False):
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
    for domain, ids in domain2dids.items():
        d = domain.lower()
        metadata = join(metadata_files_dir, f"genbank_{d}_assembly_summary.txt")
        tqdm.write(f'read {metadata}')
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
    for d in ["bacteria",'archaea']:
        metadata = join(metadata_files_dir, f"genbank_{d}_assembly_summary.txt")
        tqdm.write(f'read {metadata}')
        for row in tqdm(open(metadata)):
            if row.startswith("GC"):
                rows = row.split('\t')
                if int(rows[0]) in ids_list:
                    collect_info.append(row)
                    domain2aids[d].append(rows[0])
    return domain2aids, collect_info

# cids,cinfo = from_name2ids("Verrucomicrobia")

def main(name=None,
         odir=None,
         formats='fasta',
         ids_list=None):
    # name = "Nitrospirae;"
    # formats = 'fasta,protein-fasta'
    # odir = '/share/home-user/thliao/data/NCBI_genbank'
    formats = formats.split(',')

    odir = realpath(odir)
    if ids_list:
        domain2aids, cinfos = id2domain_to_ids(ids_list)
    else:
        domain2aids, cinfos = from_name2ids(name)

    # filter with existing files
    downloaded_aids = []
    new_domain2aids = {}
    for d, aids in domain2aids.items():
        old_d = aids[::]
        curr_dir = join(db_dir,'genbank',d)
        if formats == 'fasta':
            # check whether other kinds of files have been downloaded
            sub_aids = [_ for _ in aids
                    if not exists(join(curr_dir,_))]
            new_domain2aids[d] = sub_aids
        downloaded_aids.extend(new_domain2aids[d])
        print(f"domain: {d}, original number of ids: {len(old_d)}, now ids: {len(new_domain2aids[d])} ")
    
    
    ngd.download(**{"assembly_accessions": downloaded_aids,
                    "dry_run": False,
                    "section": "genbank",
                    "output": db_dir, # all genomes were downloaded to db_dir
                    "file_formats": formats})
    with open(join(odir,'metadata.csv'),'w') as f1:
        f1.write('\n'.join(cinfos))

@click.command()
@click.option("-n", "name", help="input the phylum name or other. use ; to separate multiple ")
@click.option("-F", "formats", help='Which formats to download (default: %(default)s).'
                                    'A comma-separated list of formats is also possible. For example: "fasta,assembly-report". '
                                    'Choose from: {choices}'.format(choices=NgdConfig.get_choices('file_formats')),
              default=NgdConfig.get_default('file_formats'))
@click.option("-o", "odir", help=f"Create output hierarchy in specified folder (default: {NgdConfig.get_default('output')}s)",
              default=NgdConfig.get_default('output'))
def cli(name, odir, formats):
    main(name, odir, formats)


if __name__ == '__main__':
    cli()
