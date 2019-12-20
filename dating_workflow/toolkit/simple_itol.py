import os
import sys

from ete3 import NCBITaxa, Tree

from api_tools.itol_func import *

ncbi = NCBITaxa()


def convert_genome_ID(genome_ID):
    # for GCA_900078535.2
    # it will return
    return genome_ID.split('_')[-1].replace('.', 'v')


def convert_genome_ID_rev(genome_ID):
    # for 900078535v2
    # it will return
    return genome_ID.replace('v', '.')


if len(sys.argv) != 3:
    raise Exception()
tree = sys.argv[1]
metadata = sys.argv[2]

gids = Tree(tree).get_leaf_names()

# leafid2gids = {_:convert_genome_ID_rev(_)
#                for _ in leafids}

# gids = set(leafid2gids.values())

gid2taxons_str = {}
gid2name = {}
gid2taxid = {}
for row in open(metadata):
    if not row.startswith("assembly_accession"):
        row = row.split('\t')
        if row[0] in gids:
            gid2name[row[0]] = row[7] + ' ' + row[9]
            gid2taxid[row[0]] = row[5]

request_taxon = 'phylum'
gid2taxon = {}
for g, tid in gid2taxid.items():
    lineage = ncbi.get_lineage(tid)
    ranks = ncbi.get_rank(lineage)
    names = ncbi.get_taxid_translator(lineage)
    r2n = {ranks[_id]: names[_id] for _id in lineage}
    gid2taxon[g] = r2n.get('phylum', '')
    if r2n.get('phylum', '') == 'Proteobacteria' and r2n.get('class', ''):
        gid2taxon[g] = r2n['class']
    gid2taxons_str[g] = ';'.join([names[_] for _ in lineage[2:]])
gid2taxon = {k: v for k, v in gid2taxon.items() if v}

text = to_label(gid2name)
os.makedirs('./itol_txt', exist_ok=1)
with open('./itol_txt/names.txt', 'w') as f1:
    f1.write(text)

text = to_label(gid2taxons_str)
os.makedirs('./itol_txt', exist_ok=1)
with open('./itol_txt/taxons_names.txt', 'w') as f1:
    f1.write(text)

text = to_label({k: k for k in gid2name})
with open('./itol_txt/reset_names.txt', 'w') as f1:
    f1.write(text)