"""
tmp file for nitrification project
"""
import os
from os.path import *
from ete3 import Tree
import re
import click
from dating_workflow.step_script import process_path
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

indir = '/home-user/thliao/template_txt/'
dataset_symbol_template = join(indir,'dataset_symbols_template.txt')



def main(intree_ori,mcmc_out_tree,output_dating_result_tree,itol_annotate,root_with):
    tree2 = Tree(intree_ori)
    tree2.set_outgroup(tree2.get_common_ancestor(root_with))
    
    mcmc_out_tree_text = open(mcmc_out_tree)
    for row in mcmc_out_tree_text:
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
    with open(output_dating_result_tree,'w') as f1:
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
    template = open('/home-user/thliao/template_txt/dataset_text_template.txt').read()
    with open(join(itol_annotate,'dating_tree_ages.txt'),'w') as f1:
        f1.write(template+'\n'+'\n'.join(raw_text))


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
    with open(join(itol_annotate,'dating_tree_bootstrap.txt'),'w') as f1:
        f1.write(template_text + annotate_text)
    
    

@click.command()
@click.option('-i','intree_ori')
@click.option('-i2','mcmc_out_tree')
@click.option('-o','output_dating_result_tree')
@click.option('-od','itol_annotate',default='./itol_txt')
@click.option('-r','root_with',help='multiple genes could use comma to separate them. LCA would be searched and taken as outgroup')
def cli(intree_ori,mcmc_out_tree,output_dating_result_tree,itol_annotate,root_with,force):
    output_dating_result_tree = process_path(output_dating_result_tree)
    itol_annotate = process_path(itol_annotate)
    if ',' in root_with:
        root_with = [_.strip() for _ in root_with.split(',')]
    else:
        root_with = [root_with.strip()]
    
    if not os.path.exists(dirname(output_dating_result_tree)):
        os.makedirs(output_dating_result_tree)
    main(intree_ori,mcmc_out_tree,output_dating_result_tree,itol_annotate,root_with)
    

# intree_ori = ''
# root_with = ''
# mcmc_out_tree = ''
# output_dating_result_tree = ''
# itol_annotate = './itol_txt'

if __name__ == "__main__":
    cli()