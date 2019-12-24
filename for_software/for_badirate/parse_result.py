import pandas as pd
from ete3 import Tree
from tqdm import tqdm

infile = "./pruned_over90Comple_lower30Contain.result"

l = '##NODES-INTERNAL_ID ASSOCIATION'
l2 = "#Branch_Group\tBirth\tDeath\tInnovation"
l3 = "##Family Turnover Rates"
l4 = "##Ancestral Family Size"
# parse, collect informations
rows = [_.strip('\n').strip('\t') for _ in open(infile).readlines()]
tree_formatted = rows[rows.index(l) + 1]
whole_tree = Tree(tree_formatted, format=1)
for n in whole_tree:
    if n.is_leaf():
        n.name = n.name.rpartition('_')[0]

fsize_idx = rows.index(l4)
gene2tree = {}
for row in rows[fsize_idx + 2:]:
    if not row:
        break
    cur_rows = row.split('\t')
    gene2tree[cur_rows[0]] = cur_rows[1]


def tree2tab(t):
    n2v = {}
    for n in Tree(t,format=1).traverse():
        if not n.name:
            leaf_names = [_.rpartition('_')[0] for _ in n.get_leaf_names()]
            n.name = whole_tree.get_common_ancestor(leaf_names).name
        if '_' in n.name:
            n.name,_,num_g = n.name.rpartition('_')
        else:
            num_g = int(n.name)
            leaf_names = [_.rpartition('_')[0] for _ in n.get_leaf_names()]
            n.name = whole_tree.get_common_ancestor(leaf_names).name
        n2v[n.name] = int(num_g)
    return n2v


# transforming into the collected
gene2n2v = {}
for gene, text_tree in tqdm(gene2tree.items()):
    if gene != 'Total Ancestral Size':
        n2v = tree2tab(text_tree)
        gene2n2v[gene] = n2v
result_df = pd.DataFrame.from_dict(gene2n2v, orient='index')
result_df.to_csv('./parsed_result.tab',sep='\t',index=1,index_label='K number')

from bin.transform.classify_kos import ko_classified_br,get_ko_infos,get_br_info
# calculating the number of genes transformed
group_dict = {"anammox":"GCA_003551305.1|GCA_004351875.1"  
              
              }
extra_tree = '../trees/iqtree/over20p_bac120.formatted.newick'


for gnum,input_str in group_dict.items():
    # input_str =   # use ';' to separate the represented IDs from two clade need to be compared.
    assert '|' in input_str
    leaf_ids,right_ids = input_str.split('|')
    all_ids = [leaf for _id in input_str.split('|') for leaf in _id.split(',')]
    LCA = whole_tree.get_common_ancestor(all_ids)
    descendent_node = whole_tree.get_common_ancestor(right_ids.split(','))
    
    LCA = LCA.name
    descendent_node = descendent_node.name
    sub_tab = result_df.loc[:, [LCA, descendent_node]]
    new_df = pd.DataFrame(index=sub_tab.index, columns=['Gain', 'Loss', 'gain/loss'])

    change = sub_tab[LCA] - sub_tab[descendent_node]
    new_df.Gain = [0 if _ <= 0 else _ for _ in change]
    new_df.Loss = [0 if _ >= 0 else _ for _ in change]
    new_df.loc[:, 'gain/loss'] = change
    new_df = new_df.loc[new_df.sum(1)!=0,:]
    
    br_kos = ko_classified_br(new_df.index)
    info_kos = get_ko_infos(new_df.index)
    infos = get_br_info(br_kos)
    {ko:kd.update({'des':info_kos.get(ko,'')}) for ko,kd in infos.items()}
    for ko in info_kos:
        if ko not in infos:
            infos[ko] = {'des':info_kos[ko]}
    ko_info_df = pd.DataFrame.from_dict(infos,orient='index')
    new_df = new_df.join(ko_info_df)
    new_df = new_df.sort_values(['gain/loss','top','A','B','C','des'])
    new_df.to_csv(f'./target_tran/group_{gnum}.csv',index_label='K number')
