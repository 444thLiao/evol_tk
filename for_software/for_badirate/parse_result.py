import pandas as pd
from ete3 import Tree
from tqdm import tqdm

infile = "./pruned_over90Comple_lower30Contain.result"

extra_tree = '../trees/iqtree/over20p_bac120.formatted.newick'
extra_t = Tree(extra_tree,format=3)

l = '##NODES-INTERNAL_ID ASSOCIATION'
l2 = "#Branch_Group\tBirth\tDeath\tInnovation"
l3 = "##Family Turnover Rates"
l4 = "##Ancestral Family Size"
# parse, collect informations
rows = [_.strip('\n').strip('\t') for _ in open(infile).readlines()]
tree_formatted = rows[rows.index(l) + 1]
whole_tree = Tree(tree_formatted, format=1)
for n in whole_tree.traverse():
    if n.is_leaf():
        n.name = n.name.rpartition('_')[0]

mapping_dict = {}
for n in whole_tree.traverse():
    if not n.is_leaf():
        all_child = n.get_leaf_names()
        LCA = extra_t.get_common_ancestor(all_child)
        mapping_dict[n.name] = LCA.name

def change_name(n):
    if n in mapping_dict:
        return mapping_dict[n]
    else:
        return n
         
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
result_df.columns = list(map(change_name,result_df.columns))
result_df.to_csv('./parsed_result.tab',sep='\t',index=1,index_label='K number')

from bin.transform.classify_kos import ko_classified_br,get_ko_infos,get_br_info
# calculating the number of genes transformed
group_dict = {"anammox":"GCA_003551305.1|GCA_001828565.1,GCA_000987375.1"  
              
              }

for gnum,input_str in group_dict.items():
    # input_str =   # use ';' to separate the represented IDs from two clade need to be compared.
    assert '|' in input_str
    leaf_ids,right_ids = input_str.split('|')
    all_ids = [leaf for _id in input_str.split('|') for leaf in _id.split(',')]
    LCA = whole_tree.get_common_ancestor(all_ids)
    descendent_node = whole_tree.get_common_ancestor(right_ids.split(','))
    
    LCA = change_name(LCA.name)
    descendent_node = change_name(descendent_node.name)
    print(f'From {LCA} to {descendent_node}')
    sub_tab = result_df.loc[:, [LCA, descendent_node]]
    
    new_df = pd.DataFrame(index=sub_tab.index, columns=['ori', 'aft', 'gain/loss'])

    change =  sub_tab[descendent_node] - sub_tab[LCA]
    new_df.ori = sub_tab[LCA]
    new_df.aft = sub_tab[descendent_node]
    new_df.loc[:, 'gain/loss'] = change
    new_df = new_df.loc[new_df.loc[:,'gain/loss']!=0,:]
    
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
    
    
    ogroup = ['GCA_001746835.1',
 'GCA_000019965.1',
 'GCA_000019665.1',
 'GCA_001318295.1',
 'GCA_001613545.1',
 'GCA_000172155.1',
 'GCA_000020225.1',
 'GCA_900097105.1',]
    ## group_test
    from scipy.stats import fisher_exact
    from statsmodels.stats.multitest import multipletests
    all_leafs = whole_tree.get_leaf_names()
    before_group = all_leafs[all_leafs.index('GCA_003551305.1'):]
    before_group = list(set(before_group).difference(set(ogroup)))
    des_group = all_leafs[all_leafs.index(all_ids[2]):all_leafs.index(all_ids[1])]
    ko2tab = {}
    ko2odd_ratio = {}
    distinct_ko = {}
    for ko in tqdm(result_df.index):
        
        ko_b_1 = sum(result_df.loc[ko,before_group]>0)
        ko_b_0 = sum(result_df.loc[ko,before_group]==0)
        ko_d_1 = sum(result_df.loc[ko,des_group]>0)
        ko_d_0 = sum(result_df.loc[ko,des_group]==0)
        tab_num = [[ko_b_1,ko_b_0],[ko_d_1,ko_d_0]]
        ko2tab[ko] = tab_num
        oddratio, p = fisher_exact(tab_num, alternative="two-sided")
        distinct_ko[ko] = p
        ko2odd_ratio[ko] = oddratio

    fe_corrected_p = multipletests([_ for k, _ in distinct_ko.items()],
                                method='fdr_bh')
    fe_corrected_ko2p = dict(zip([k for k, _ in distinct_ko.items()],
                                fe_corrected_p[1]))
    sig_ko_list = {k: v for k, v in fe_corrected_ko2p.items() if v <= 0.05}
    len(new_df.index.intersection(set(sig_ko_list)))
    
    new_df.loc[new_df.index.intersection(set(sig_ko_list)),'significance'] = 'True'
    new_df.to_csv(f'./target_tran/group_{gnum}.csv',index_label='K number')

from api_tools.itol_func import to_binary_shape

sig_ko = list(new_df.index.intersection(set(sig_ko_list)))[-30:]
all_gids = whole_tree.get_leaf_names()
id2ko = {}
for gid in all_gids:
    
    remained_kos = set(result_df.index[result_df.loc[:,gid] > 0]).intersection(set(sig_ko))
    id2ko[gid] = list(remained_kos)
text = to_binary_shape(id2ko,)
with open('./target_tran/test_a_p.txt','w') as f1:
    f1.write(text)