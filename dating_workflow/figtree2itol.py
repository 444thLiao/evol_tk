"""
tmp file for nitrification project
"""
import os
from os.path import *
from ete3 import Tree
import re
#   Tree with NHX style metadata:

#    (A:0.1,(B:0.2,(C:0.2,D:0.3):0.4[&&NHX:conf=0.01:name=NODE1]):0.5);

#     A, B, C   : leaf names
#     internal node will have the ID NODE1
#     metadata value 'conf' will be available for visualization

def sub_for(m):
    t = m.string[m.start(0):m.end(0)]
    t = t.replace(',','-')
    #t = t.replace('.','v')
    t = t.replace(' ','')
    t = t.replace('[&95%HPD={',"")
    t = t.replace('}]',"")
    #t = t.replace(',','_')
    #print(t.string)
    return t

mcmc_out_tre = './FigTree.newick'
a = open(mcmc_out_tre)
for row in a:
    if row.strip().startswith('UTREE 1 ='):
        t = row.split('UTREE 1 =')[1].strip('\n')
        t = re.sub('\[&95%HPD=.*?\]',sub_for,t)
        t = t.replace(' ','')
        tree = Tree(t,format=1)
        
count = 0 
for n in tree.traverse():
    if not n.is_leaf():
        dates = n.name
        n.name = 'I%s' % count
        n.add_features(ages=dates,)
        all_leafs = n.get_leaf_names()
        nin2 = tree2.get_common_ancestor(all_leafs)
        n.add_features(support=nin2.support)
        count +=1
# tree.features.remove('support')
text = tree.write(format=3)
text = text.replace(')1:','):')
with open('./manual_itol.nexus','w') as f1:
    f1.write(text)
    
    
raw_text = []
for n in tree.traverse():
    if not n.is_leaf():
        raw_text.append("\t".join([n.name,
                                   n.ages,
                                   '1',
                                   	'#FF0000',
                                    'bold',
                                    '1',
                                    '0'
                               ]))
with open('./itol_txt/itol_ages.txt','w') as f1:
    f1.write('\n'.join(raw_text))

    
indir = '/home-user/thliao/template_txt/'
dataset_symbol_template = join(indir,'dataset_symbols_template.txt')

root_with = '''GCA_000011385.1
GCA_000013205.1
GCA_000317065.1
GCA_000332175.1
GCA_000332215.1
GCA_000011345.1
GCA_000022045.1
GCA_000018105.1'''.split('\n')
ori_tre = "/home-user/thliao/data/nitrification_for/dating_for/ppostcluster/iqtree.treefile"
tree2 = Tree(ori_tre)
tree2.set_outgroup(tree2.get_common_ancestor(root_with))
rows = []
template_text = open(dataset_symbol_template).read()
for n in tree2.traverse():
    if not n.is_leaf():
        size = '5'
        shape = '2'
        filled = '1'
        s_v = n.support
        childrens = n.get_leaf_names()
        nid = tree.get_common_ancestor(childrens)
        
        if int(s_v) >= 95:
            color = '#000000'
        elif int(s_v) >= 85 and int(s_v) < 95:
            color = '#777777'
        elif int(s_v) >= 65 and int(s_v) < 85:
            color = '#eeeeee'
        else:
            color = '#FFFFFF'
            
        if color:
            row = '\t'.join([nid.name,shape,size,color,filled,'0',''])
            rows.append(row)

    annotate_text = '\n'.join(rows)
    template_text = template_text.format(dataset_label='bootstrap',
                                         legend_text='',
                                         maximum_size=10)
with open('./itol_txt/bootstrap.txt','w') as f1:
    f1.write(template_text + annotate_text)
