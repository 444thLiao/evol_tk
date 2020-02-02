"""
This script is considerated the backbone gene presence situtations

"""



# annotate 27 genes
from Bio import SeqIO
from dating_workflow.step_script.extract_cog25 import parse_annotation
from api_tools.itol_func import *
# rrna_dir = './rrna'
# gid2genes = {k: [_k for _k, _v in v.items() if _v] for k, v in _subgenome2cdd.items()}

# for record in SeqIO.parse(join(rrna_dir, '16S.fasta'), format='fasta'):
#     gname = 'GCA_' + convert_genome_ID_rev(record.id.split('_')[0])
#     if gname in gid2genes:
#         gid2genes[gname].append('16S')
# for record in SeqIO.parse(join(rrna_dir, '23S.fasta'), format='fasta'):
#     gname = 'GCA_' + convert_genome_ID_rev(record.id.split('_')[0])
#     if gname in gid2genes:
#         gid2genes[gname].append('23S')

in_annotations = './cog25_annotate'
evalue = 1e-20
gid2genes = parse_annotation(in_annotations, top_hit=True,evalue=evalue)
all_genes = set([_ for vl in gid2genes.values() for _ in vl])
text = to_binary_shape(gid2genes,
                       {g: {'color': '#007acc'} for g in all_genes})

with open('./itol_txt/25genes.txt', 'w') as f1:
    f1.write(text)

from ete3 import Tree
import pandas as pd
tree = './trees/final/198g_merged.newick'
t = Tree(tree,format=3)
all_ids = t.get_leaf_names()

id2gene = {genome: {k:1 if v else 0 for k,v in _v.items() } 
           for genome,_v in gid2genes.items() }
info_df = pd.DataFrame.from_dict(id2gene,orient='index')
info_df = info_df.reindex(all_ids)
info_df.to_excel('./itol_txt/25genes.xlsx',index=1)

# ids = ['CDD:223275','CDD:223172','CDD:223280']
# gids = info_df.index[(info_df.loc[:,ids] !=1).any(1)]

# text = to_binary_shape({g:'1' for g in gids},
#                        {'1': {'color': '#007acc'}})
# with open('./itol_txt/tmp.txt', 'w') as f1:
#     f1.write(text)
ids = open('./dating_for_190g.list').read().split('\n')
text = to_binary_shape({k:['remained'] for k in ids})
with open('./itol_txt/test_remained.txt','w') as f1:
   f1.write(text)
   
id2gene = {genome: {k:1 if v else 0 for k,v in _v.items() } 
           for genome,_v in gid2genes.items() if genome in ids}
info_df = pd.DataFrame.from_dict(id2gene,orient='index')
info_df = info_df.reindex([_ for _ in all_ids if _ in ids])
info_df.to_excel('./itol_txt/25genes.xlsx',index=1)
