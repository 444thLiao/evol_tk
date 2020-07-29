"""
For implement some missing feature(but i need) for ncbi-genome-download
"""

from os.path import *

import click
import ncbi_genome_download as ngd
from ete3 import NCBITaxa
from ncbi_genome_download import NgdConfig
from tqdm import tqdm


def from_name2ids(phylum_name):
    phylum_names = [_ for _ in phylum_name.split(';') if _]
    # phylum_name = "Nitrospirae;"
    # phylum_tid = "40117"
    ncbi = NCBITaxa()

    tid = ncbi.get_name_translator(phylum_names)

    for _ in phylum_names:
        if not tid.get(_):
            print(f" '{_}'' not found. please check the name")
    tids = [tid.get(_, [None])[0]
            for _ in phylum_names
            if tid.get(_)]

    descend_ids = []
    for tid in tids:
        _descend_ids = ncbi.get_descendant_taxa(tid,
                                                intermediate_nodes=True)
        descend_ids += _descend_ids
    print(f"in total, {len(descend_ids)} taxids were found. ")

    metadata = expanduser("~/.cache/ncbi-genome-download/genbank_bacteria_assembly_summary.txt")
    collect_ids = []
    collect_info = []
    for row in tqdm(open(metadata)):
        if row.startswith("GC"):
            rows = row.split('\t')
            if int(rows[5]) in descend_ids:
                collect_info.append(row)
                collect_ids.append(rows[0])
    return collect_ids, collect_info
# cids,cinfo = from_name2ids("Verrucomicrobia")

def main(name, odir, formats, ids_list):
    # name = "Nitrospirae;"
    # formats = 'fasta,protein-fasta'
    # odir = '/share/home-user/thliao/data/NCBI_genbank'
    formats = formats.split(',')

    odir = realpath(odir)
    cids, cinfos = from_name2ids(name)
    ngd.download(**{"assembly_accessions": cids,
                    "dry_run": False,
                    "section": "genbank",
                    "output": odir,
                    "file_formats": formats})


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
