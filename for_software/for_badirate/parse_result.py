import pandas as pd
from ete3 import Tree

infile = "./test.out"

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

# calculating the number of genes transformed
input_str = 'dyak;dmel'  # use ';' to separate the represented IDs from two clade need to be compared.
LCA = whole_tree.get_common_ancestor(input_str.split(';'))
LCA_anc = LCA.get_ancestors()[0]

LCA = LCA.name
LCA_anc = LCA_anc.name
sub_tab = result_df.loc[:, [LCA, LCA_anc]]
new_df = pd.DataFrame(index=sub_tab.index, columns=['Gain', 'Loss', 'gain/loss'])

change = sub_tab[LCA] - sub_tab[LCA_anc]
new_df.Gain = [0 if _ <= 0 else _ for _ in change]
new_df.Loss = [0 if _ >= 0 else _ for _ in change]
new_df.loc[:, 'gain/loss'] = change
new_df.to_csv()
