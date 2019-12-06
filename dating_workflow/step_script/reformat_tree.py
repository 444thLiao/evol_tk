"""
reformat the result tree from formatted into tree without branch and with calibrations.
"""
from ete3 import Tree

tree_file = '../trees/iqtree/170g_concat.formatted.newick'
calibration_txt = './calibration.txt'
# root_at = 'GCA_000011385.1,GCA_000009725.1'


calibration_dict = {}
for row in open(calibration_txt):
    if row and not row.startswith('#'):
        LCA,time,remained = row.split('\t')
        calibration_dict[LCA] = time
        

t = Tree(tree_file,format=3)
# new_root = t.get_common_ancestor([l 
#                                   for l in t.get_leaves() 
#                                   if l.name in root_at])
# t.set_outgroup(new_root)

for LCA,time in calibration_dict.items():
    names = LCA.split('|')
    LCA_node = t.get_common_ancestor([l
                                      for l in t.get_leaves()
                                      if l.name in names])
    LCA_node.name = time

for n in t.traverse():
    if not n.is_leaf():
        if 'I' in n.name and 'S' in n.name:
            n.name = ''
final_tree = t.write(format=8,format_root_node=True).replace('NoName','')

with open('./170g_bac120_3cal.newick','w') as f1:
    f1.write(f'{len(t.get_leaves())}\n'+final_tree)
    
    
from api_tools.itol_func import *

size = '12'
shape = '2'
filled = '1'
color = '#ff9900'

calibration_txt = './calibration.txt'
calibration_dict = {}

rows_str = []
rows = []
for row in open(calibration_txt):
    if row and not row.startswith('#'):
        LCA,time,remained = row.split('\t')
        
        row = '\t'.join([LCA,shape,size,color,filled,'1',time])
        rows.append(row)
        row = '\t'.join([LCA,time,'1','#FF0000','bold','2','0'])
        rows_str.append(row)
template_text = open(dataset_symbol_template).read()
annotate_text = '\n'.join(rows)
template_text = template_text.format(dataset_label='calibration',
                                        legend_text='',
                                        maximum_size=size)
with open('../itol_txt/itol_calibration.txt','w') as f1:
    f1.write(template_text+'\n'+annotate_text)

template_text = open('/home-user/thliao/template_txt/dataset_text_template.txt').read()
annotate_text = '\n'.join(rows_str)
with open('../itol_txt/itol_calibration_str.txt','w') as f1:
    f1.write(template_text+'\n'+annotate_text)

text = to_node_symbol(tree_file)
with open('../itol_txt/bootstraps.txt','w') as f1:
    f1.write(text)