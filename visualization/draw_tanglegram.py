# from my_toolkit

"python3 draw_tanglegram.py -newick1 ./all_1469_new.newick -newick2 ./nxrA.newick -cf1 ./gene_annotation.txt -cf2 ./phylum_annotate.txt -length 'max' -sep '_' -extra_set 'rename' "
import sys

from .tanglegram import *
from bin.format_newick import sort_tree
from os.path import dirname, join, exists
from ete3 import Tree
import io



# raw gene 2 species tree
example_dir = r"D:\Desktop\OneDrive - The Chinese University of Hong Kong\luo lab\项目\AOB\whole_tree\nxrA"
desktop_NOB_tmp = "/d/Desktop/NOB_HGT"
gene_tree = join(example_dir, 'nxrA.newick')
species_tree = join(example_dir,'..', 'all','all_1469_new.newick')
species_tree = Tree(species_tree,format=3)
for l in species_tree.get_leaves():
    l.name = l.name.split('_')[-1].replace('.','v')



gene_tree_colors = join(example_dir, 'gene_annotation.txt')
species_tree_colors = join(example_dir,'..', 'all', 'phylum_annotate.txt')

Angst_reconciles = sort_tree(Tree("D:/Desktop/NOB_HGT/angst/AnGST.newick",format=3))

fig = main(gene_tree,species_tree,
           gene_tree_colors, species_tree_colors,
           l_legnth='max', sep='_', extra_set='rename')
c2name = parse_color_scheme_files(gene_tree_colors,get_raw_name=True)
middle_lines = [_ for _ in fig.data if _['xaxis']=='x2']
for t in middle_lines:
    t.name = c2name.get(t['line']['color'], "no reference")
    
fig.layout.width = 1500
fig.layout.height = 4000
fig.layout.yaxis.zeroline = False
fig.layout.yaxis2.zeroline = False
fig.layout.yaxis3.zeroline = False
fig.layout.template = 'xgridoff'
fig.write_html(join(example_dir,'tanglegram', 'gene2species.html'))

# Angst(reconcile gene) 2 species tree

fig = main(Angst_reconciles,species_tree,
           gene_tree_colors, species_tree_colors,
           l_legnth='max', sep='_', extra_set='rename')
middle_lines = [_ for _ in fig.data if _['xaxis']=='x2']
for t in middle_lines:
    t.name = c2name.get(t['line']['color'], "no reference")
fig.layout.width = 1500
fig.layout.height = 4000
fig.layout.yaxis.zeroline = False
fig.layout.yaxis2.zeroline = False
fig.layout.yaxis3.zeroline = False
fig.layout.template = 'xgridoff'
fig.write_html(join(example_dir,'tanglegram', 'Angst2species.html'))

# Angst(reconcile gene) 2 raw gene

fig = main(Angst_reconciles,gene_tree,
           gene_tree_colors, gene_tree_colors,
           l_legnth='max', sep='_', extra_set=False,identical=True)
middle_lines = [_ for _ in fig.data if _['xaxis']=='x2']
for t in middle_lines:
    t.name = c2name.get(t['line']['color'], "no reference")
fig.layout.width = 1500
fig.layout.height = 4000
fig.write_html(join(example_dir,'tanglegram', 'Angst2gene.html'))


# Angst(reconcile gene) 2 truncated genome tree
truncated_species_tree = r"D:\Desktop\NOB_used_first_pruned.newick"
truncated_species_tree = Tree(truncated_species_tree,format=3)
for l in truncated_species_tree.get_leaves():
    l.name = l.name.split('_')[-1].replace('.','v')

fig = main(Angst_reconciles,truncated_species_tree,
           gene_tree_colors, species_tree_colors,
           l_legnth='max', sep='_', extra_set='rename')
middle_lines = [_ for _ in fig.data if _['xaxis']=='x2']
for t in middle_lines:
    t.name = c2name.get(t['line']['color'], "no reference")
fig.layout.width = 1500
fig.layout.height = 4000
fig.write_html(r"D:\Desktop\test.html")
fig.write_html(join(example_dir,'tanglegram', 'Angst2gene.html'))



# dated  (ecceTERA)
species_tree_159g = './NOB_HGT/159g_nucl_template/159g_dating.newick'
species_tree_159g = Tree(species_tree_159g,format=3)
for l in species_tree_159g.get_leaves():
    l.name = l.name.split('_')[-1].replace('.','v')
    
ecceTERA_gene_tree = sort_tree('./NOB_HGT/ecceTERA_2/geneTree')

fig = main(ecceTERA_gene_tree,species_tree_159g,
           gene_tree_colors, species_tree_colors,
           l_legnth='max', sep='_', extra_set='rename')
middle_lines = [_ for _ in fig.data if _['xaxis']=='x2']
for t in middle_lines:
    t.name = c2name.get(t['line']['color'], "no reference")
fig.layout.yaxis.zeroline = False
fig.layout.yaxis2.zeroline = False
fig.layout.yaxis3.zeroline = False
fig.layout.width = 1000
fig.layout.height = 3000
fig.layout.template = 'xgridoff'
fig.write_html(r"D:\Desktop\test.html")


# dated  (ecceTERA)
species_tree_159g = './NOB_HGT/159g_nucl_template/159g_dating.newick'
species_tree_159g = Tree(species_tree_159g,format=3)
for l in species_tree_159g.get_leaves():
    l.name = l.name.split('_')[-1].replace('.','v')
    
ecceTERA_gene_tree = sort_tree('./NOB_HGT/angst_2/geneTree')

fig = main(ecceTERA_gene_tree,species_tree_159g,
           gene_tree_colors, species_tree_colors,
           l_legnth='max', sep='_', extra_set='rename')
middle_lines = [_ for _ in fig.data if _['xaxis']=='x2']
for t in middle_lines:
    t.name = c2name.get(t['line']['color'], "no reference")
fig.layout.yaxis.zeroline = False
fig.layout.yaxis2.zeroline = False
fig.layout.yaxis3.zeroline = False
fig.layout.width = 1000
fig.layout.height = 3000
fig.layout.template = 'xgridoff'
fig.write_html(r"D:\Desktop\test.html")