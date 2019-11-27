from api_tools.itol_func import *
import sys
from ete3 import NCBITaxa,Tree
import os
ncbi = NCBITaxa()

def convert_genome_ID(genome_ID):
    # for GCA_900078535.2
    # it will return
    return genome_ID.split('_')[-1].replace('.','v')

def convert_genome_ID_rev(genome_ID):
    # for 900078535v2
    # it will return
    return genome_ID.replace('v','.')

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
            gid2name[row[0]] = row[7]
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


id2info = gid2taxon
id2info, info2color = get_colors_general(id2info)
text = to_color_strip(id2info, info2color, info_name='phylum')

with open('./itol_txt/phylum_annotate.txt', 'w') as f1:
    f1.write(text)

id2info = gid2taxon
id2info, info2color = get_colors_general(id2info)
text = to_color_branch(id2info, info2color, dataset_name='phylum/class',no_legend=True)

with open('./itol_txt/phylum_annotate_branch.txt', 'w') as f1:
    f1.write(text)
    
    
# annotate 27 genes
gid2genes = {k:[_k for _k,_v in v.items() if _v] for k,v in genome2cdd.items()}
text = to_binary_shape(gid2genes)
