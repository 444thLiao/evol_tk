import pandas as pd
from ete3 import Tree
from bin.other_convertor.classify_kos import *

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

def fe_text(g1, g2, ko_df):
    ko2tab = {}
    ko2odd_ratio = {}
    ko2pvalue = {}
    for ko in tqdm(set(ko_df.columns)):
        ko_g1_p = sum(ko_df.loc[g1, ko] > 0)
        ko_g1_a = sum(ko_df.loc[g1, ko] == 0)
        ko_g2_p = sum(ko_df.loc[g2, ko] > 0)
        ko_g2_a = sum(ko_df.loc[g2, ko] == 0)
        tab_num = [[ko_g1_p, ko_g1_a],
                   [ko_g2_p, ko_g2_a]]
        ko2tab[ko] = tab_num
        oddratio, p = fisher_exact(tab_num, alternative="two-sided")
        ko2pvalue[ko] = p
        ko2odd_ratio[ko] = oddratio

    fe_corrected_p = multipletests([_ for k, _ in ko2pvalue.items()],
                                   method='fdr_bh')
    fe_corrected_ko2p = dict(zip([k for k, _ in ko2pvalue.items()],
                                 fe_corrected_p[1]))
    sig_ko_list = {k: v for k, v in fe_corrected_ko2p.items() if v <= 0.05}
    return ko2tab,ko2odd_ratio,ko2pvalue,sig_ko_list


# all gene table
tab = "/mnt/home-backup/thliao/plancto/protein_annotations/hmmsearch_merged/merged_hmm_binary.tab"
all_df = pd.read_csv(tab, sep='\t', index_col=0)

# for dating tree
intree = '/mnt/home-backup/thliao/plancto/trees/final/83g_merged_sorted.newick'
t = Tree(intree, 3)

## anammox group
LCA = t.get_common_ancestor(['GCA_001828545.1', 'GCA_001828295.1'])
larger_LCA = t.get_common_ancestor(['GCA_003694635.1', 'GCA_003551565.1'])

#
g = larger_LCA.get_leaf_names()
g1 = LCA.get_leaf_names()
# target group
g2 = set(g).difference(set(g1))
# remainning group
ko2tab,ko2odd_ratio,ko2pvalue,sig_ko_list = fe_text(g1,g2,all_df)

# summarized df
br_kos = ko_classified_br(sig_ko_list)
md_kos = ko_classified_module(sig_ko_list)
md2info = get_md_infos(md_kos.keys())
info_kos = get_ko_infos(sig_ko_list)
infos = get_br_info(br_kos)
df = pd.concat(infos, axis=0)
df.loc[:, 'des'] = [info_kos.get(_, '') for _ in df.index]

new_df = pd.DataFrame(index=sig_ko_list, columns=['num_p in target', 'num_a in target',
                                                  'num_p in remained', 'num_a in remained',
                                                  ])
for ko in new_df.index:
    new_df.loc[ko,:] = ko2tab[ko][0] + ko2tab[ko][1]
new_df = new_df.reindex(index=df.index)

final_df = pd.concat([new_df, df], axis=1)
for md, kos in md_kos.items():
    final_df.loc[kos, "module ID"] = md
    final_df.loc[kos, "module name"] = md2info[md]['name']
    final_df.loc[kos, "module all ko"] = md2info[md]['all ko']

final_df = final_df.sort_values(['num_p in target', 'top', 'A', 'B', 'C', 'des'])
final_df.to_excel('./test_83g_phyletic_pattern.xlsx')



## anammox inside
LCA = t.get_common_ancestor(['GCA_004282745.1', 'GCA_008636105.1'])
larger_LCA = t.get_common_ancestor(['GCA_004282745.1', 'GCA_007618145.1'])

#
g = larger_LCA.get_leaf_names()
g1 = LCA.get_leaf_names()
# target group
g2 = set(g).difference(set(g1))
# remainning group
ko2tab,ko2odd_ratio,ko2pvalue,sig_ko_list = fe_text(g1,g2,all_df)

# summarized df
br_kos = ko_classified_br(sig_ko_list)
md_kos = ko_classified_module(sig_ko_list)
md2info = get_md_infos(md_kos.keys())
info_kos = get_ko_infos(sig_ko_list)
infos = get_br_info(br_kos)
df = pd.concat(infos, axis=0)
df.loc[:, 'des'] = [info_kos.get(_, '') for _ in df.index]

new_df = pd.DataFrame(index=sig_ko_list, columns=['num_p in target', 'num_a in target',
                                                  'num_p in remained', 'num_a in remained',
                                                  ])
for ko in new_df.index:
    new_df.loc[ko,:] = ko2tab[ko][0] + ko2tab[ko][1]
new_df = new_df.reindex(index=df.index)

final_df = pd.concat([new_df, df], axis=1)
for md, kos in md_kos.items():
    final_df.loc[kos, "module ID"] = md
    final_df.loc[kos, "module name"] = md2info[md]['name']
    final_df.loc[kos, "module all ko"] = md2info[md]['all ko']

final_df = final_df.sort_values(['num_p in target', 'top', 'A', 'B', 'C', 'des'])
final_df.to_excel('./test_83g_phyletic_pattern_insideanammox.xlsx')