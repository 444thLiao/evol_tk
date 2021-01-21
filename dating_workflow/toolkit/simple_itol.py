import os
import sys

from ete3 import NCBITaxa, Tree
from os.path import *
from api_tools.itol_func import *
from dating_workflow.step_script import process_path,convert_genome_ID,convert_genome_ID_rev
ncbi = NCBITaxa()


if len(sys.argv) != 3:
    raise Exception()



genome_list = sys.argv[1]
metadata = expanduser('~/.cache/ncbi-genome-download/genbank_bacteria_assembly_summary.txt')

gids = [_ for _ in open(genome_list).read().split('\n') if _]

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

## 25 genes
from dating_workflow.step_script.extract_cog25 import parse_annotation
in_annotations = '/mnt/home-backup/thliao/NCBI/modified_data/cog25_annotate'
genome2cdd = parse_annotation(in_annotations, top_hit=True, evalue=1e-20)
genomes_list = [_ for _ in open('./over20p_genomes_add_cyano.list').read().split('\n') if _]
gid2genes = {k: [_k for _k, _v in v.items() if _v] for k, v in genome2cdd.items() if k in genomes_list}
all_genes = set([_ for vl in gid2genes.values() for _ in vl])
text = to_binary_shape(gid2genes,
                       {g: {'color': '#007acc'} for g in all_genes})

with open('./itol_txt/cog25genes.txt', 'w') as f1:
    f1.write(text)
    
    
    
## Genome type
import pandas as pd
metadata = expanduser('~/.cache/ncbi-genome-download/genbank_bacteria_assembly_summary.txt')
t = open(metadata)
next(t)
header = next(t).strip('#\n ').split('\t')
df = pd.read_csv(metadata,sep='\t',header=None,comment='#',low_memory=False)


