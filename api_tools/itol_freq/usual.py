import os

import pandas as pd
from ete3 import NCBITaxa
from tqdm import tqdm
from api_tools.itol_func import *
from dating_workflow.step_script import convert_genome_ID_rev
import plotly.express as px

ncbi = NCBITaxa()
home = os.getenv("HOME")
metadata = f"{home}/.cache/ncbi-genome-download/genbank_bacteria_assembly_summary.txt"

gids = [_
        for _ in open('./rawdata/ids.list').read().split('\n')
        if _.startswith('GCA')]
gids = set(gids)
gid2taxons_str = {}
gid2name = {}
gid2taxid = {}
for row in tqdm(open(metadata)):
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
    gid2taxon[g] = r2n.get(request_taxon, '')
    if r2n.get('phylum', '') == 'Proteobacteria' and r2n.get('class', ''):
        gid2taxon[g] = r2n['class']
    elif r2n.get('phylum', '') == "candidate division NC10":
        gid2taxon[g] = "NC10"
    gid2taxons_str[g] = ';'.join([names[_] for _ in lineage[2:]])
gid2taxon = {k: v for k, v in gid2taxon.items() if v}

# taxonomy annotation
color_scheme = {'type': {'NOB': '#e41a1c', 'comammox': '#edc31d',
                         'AOB': '#bad5b9', 'AOA': '#358f0f'},
                'phylum/class': {'Thaumarchaeota': '#358f0f',
                                 'Nitrospirae': '#edc31d',
                                 'Gammaproteobacteria': '#78fce0',
                                 'Chloroflexi': '#e41a1c',
                                 'Betaproteobacteria': '#956cb4',
                                 'Alphaproteobacteria': '#8c613c',
                                 'Actinobacteria': '#11FF11',
                                 'Planctomycetes': '#FF66bb',
                                 "Nitrospinae": "#4285f4",
                                 "Verrucomicrobia": "#E57373",
                                 "CPR": "#74a45b",
                                 "NC10": "#7986CB",
                                 "Armatimonadetes": "#C86437",
                                 }}


def get_colors_general(ID2infos, now_info2style={}):
    _c = color_scheme.copy()
    _c = _c['phylum/class']
    colors = px.colors.qualitative.Dark24 + px.colors.qualitative.Light24
    remained_colors = [c for c in colors if c not in now_info2style.values()]
    info2style = {}
    for v in set(ID2infos.values()):
        if v in _c:
            info2style.update({v: _c[v]})
        else:
            one_color = remained_colors.pop(0)
            info2style.update({v: one_color})
    return ID2infos, info2style


odir = './itol_txt'
os.makedirs('./itol_txt', exist_ok=True)
id2info, info2color = get_colors_general({k: v
                                          for k, v in gid2taxon.items()
                                          if v in color_scheme["phylum/class"]})
id2info = {k: v
               for k, v in id2info.items()
               if v in color_scheme['phylum/class']}
sub_info = {k: v
            for k, v in color_scheme["phylum/class"].items()
            if k in id2info.values()}
text = to_color_strip(id2info, sub_info, info_name='phylum')
with open(join(odir, 'phylum_annotate.txt'), 'w') as f1:
    f1.write(text)

text = to_label(gid2name)
with open('./itol_txt/names.txt', 'w') as f1:
    f1.write(text)

text = to_label(gid2taxons_str)
with open('./itol_txt/taxons_names.txt', 'w') as f1:
    f1.write(text)

text = to_label({k: k for k in gid2name})
with open('./itol_txt/reset_names.txt', 'w') as f1:
    f1.write(text)

id2info = gid2taxon
id2info, info2color = get_colors_general(id2info)
text = to_color_branch(id2info, info2color, dataset_name='phylum/class', no_legend=True)

with open('./itol_txt/phylum_annotate_branch.txt', 'w') as f1:
    f1.write(text)

# annotate 27 genes
from Bio import SeqIO

rrna_dir = './rrna'
gid2genes = {k: [_k for _k, _v in v.items() if _v] for k, v in _subgenome2cdd.items()}

for record in SeqIO.parse(join(rrna_dir, '16S.fasta'), format='fasta'):
    gname = 'GCA_' + convert_genome_ID_rev(record.id.split('_')[0])
    if gname in gid2genes:
        gid2genes[gname].append('16S')
for record in SeqIO.parse(join(rrna_dir, '23S.fasta'), format='fasta'):
    gname = 'GCA_' + convert_genome_ID_rev(record.id.split('_')[0])
    if gname in gid2genes:
        gid2genes[gname].append('23S')

all_genes = set([_ for vl in gid2genes.values() for _ in vl])
text = to_binary_shape(gid2genes,
                       {g: {'color': '#007acc'} for g in all_genes})

with open('./itol_txt/27genes.txt', 'w') as f1:
    f1.write(text)
# annotate cog25
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
    
# annotate completeness

sub_ids = [_ for _ in open('./over20p_genomes.list').read().split('\n') if _]
completeness_df = pd.read_csv('./checkM_result_phylum/merged_checkM.tab', sep='\t', index_col=0)
sub_df = completeness_df.reindex(sub_ids)
sub_df = sub_df.loc[~sub_df.isna().all(1)]
id2val = dict(zip(sub_ids, sub_df.loc[:, 'Completeness']))
text = color_gradient(id2val,
                      max_val=max(sub_df.loc[:, 'Completeness']),
                      min_val=min(sub_df.loc[:, 'Completeness']),
                      mid_val=50)
with open('./itol_txt/over20p_completeness.txt','w') as f1:
    f1.write(text)
text = get_text_anno(id2val,extra_replace={"#HEIGHT_FACTOR,1":"HEIGHT_FACTOR\t1.5",
                                           "#HORIZONTAL_GRID,1":"HORIZONTAL_GRID\t0",
                                           "#VERTICAL_GRID,1":"VERTICAL_GRID\t0",
                                           })
with open('./itol_txt/over20p_completeness_text.txt','w') as f1:
    f1.write(text)


# bootstrap
intree = ''
tree = Tree(intree,format=3)
text = to_node_symbol(tree)
with open('./itol_txt/bp_rrna.txt','w') as f1:
    f1.write(text)

