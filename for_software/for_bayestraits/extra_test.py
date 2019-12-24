"""
For testing one node which is LCA of given leaf ids.
"""

# compare LCA of both left and righ to the nearest ancestor of right
from for_software.for_bayestraits.toolkit.get_result import get_result, summaized_r, summaized_rate
import multiprocessing as mp
from for_software.for_bayestraits.toolkit.construct_kit import nw2nexus, get_tags
from tqdm import tqdm
from ete3 import Tree
import click
from subprocess import check_call
from os.path import *
import os
target_group = "GCA_002007645.1|GCA_001825665.1"  # M2N
"GCA_003576905.1|GCA_008363155.1"   # M2N
"GCA_003576845.1|GCA_002277955.1"  # M2N
"GCA_002709045.1|GCA_007751035.1"  # M2N
"GCA_007859755.1|GCA_003335505.1,GCA_002967715.1"  # M2N
"GCA_007694185.1|GCA_000186345.1,GCA_005792915.1"  # M2N
"GCA_005239935.1|GCA_007693455.1"  # root

"GCA_003551305.1|GCA_004351875.1"  # non-anammox 2 anammox
"GCA_004376375.1"  # alternative more deeper non-anammox


# intree = './trees/iqtree/over20p_bac120.formatted.newick'
# inmetadata = './bayesTraits_test/m2nm.txt'
# odir = './bayestraits_habitat'

bt_exe = expanduser("~/software/BayesTraitsV3.0.2-Linux/BayesTraitsV3")


def run(cmd):
    check_call(cmd + ' >/dev/null', shell=True)


intree = './over20p_bac120.formatted.newick'
metadata = './metadata.txt'
odir = './sig_test'

ori_f = join('./', 'complex_m', 'bst_complex.Stones.txt')

# intree = '/home-user/thliao/data/plancto/trees/iqtree/over20p_bac120.formatted.newick'
# t = Tree(intree, format=3)
# for idx, node_str in enumerate(groups_list):
#     leaf_ids, right_ids = node_str.split('|')
#     all_ids = [leaf for _id in node_str.split('|') for leaf in _id.split(',')]
#     LCA = t.get_common_ancestor(all_ids)
#     print(LCA.name)
groups_dict = {"I104": "GCA_002007645.1|GCA_001825665.1",  # M2N
               "I102": "GCA_003576905.1|GCA_008363155.1",  # M2N
               "I93": "GCA_003576845.1|GCA_002277955.1",  # M2N
               "I438": "GCA_002709045.1|GCA_007751035.1",  # M2N
               "I667": "GCA_007859755.1|GCA_003335505.1,GCA_002967715.1",  # M2N
               "I41": "GCA_007694185.1|GCA_000186345.1,GCA_005792915.1", # M2N
               "I1": "GCA_005239935.1|GCA_007693455.1",  # root
               "I3":"GCA_002050035.1|GCA_002746535.1",     
               }

groups_dict = {"I104": "GCA_002007645.1|GCA_001825665.1",  # M2N
               "I102": "GCA_003576905.1|GCA_008363155.1",  # M2N
               "I93": "GCA_003576845.1|GCA_002277955.1",  # M2N
               "I438": "GCA_002709045.1|GCA_007751035.1",  # M2N
               "I667": "GCA_007859755.1|GCA_003335505.1,GCA_002967715.1",  # M2N
               "I41": "GCA_007694185.1|GCA_000186345.1,GCA_005792915.1", # M2N
               "I1": "GCA_005239935.1|GCA_007693455.1",  # root
               "I3":"GCA_002050035.1|GCA_002746535.1",       
               "I141":"GCA_001824635.1|GCA_001304275.1",
               "I78":"GCA_001999965.1|GCA_001304275.1",
               
               }

Force_state1 = 'M'
Force_state2 = 'N'

params = []
for name, node_str in groups_dict.items():
    leaf_ids, right_ids = node_str.split('|')
    all_ids = [leaf for _id in node_str.split('|') for leaf in _id.split(',')]

    complex_model = ["1",
                     "2",
                     "PriorAll exp 10",
                     "Stones 100 1000",
                     "AddTag FNode " + ' '.join(all_ids),
                     f"Fossil Node01 FNode {Force_state1}"
                     ]

    ofile1 = join(odir, f'test_{name}_params_{Force_state1}.txt')
    with open(ofile1, 'w') as f1:
        f1.write('\n'.join(complex_model +
                           [f"LF {ofile1.replace('_params.txt','')}",
                            'Run']))
    complex_model = ["1",
                     "2",
                     "PriorAll exp 10",
                     "Stones 100 1000",
                     "AddTag FNode " + ' '.join(all_ids),
                     f"Fossil Node01 FNode {Force_state2}"
                     ]

    ofile2 = join(odir, f'test_{name}_params_{Force_state2}.txt')
    with open(ofile2, 'w') as f1:
        f1.write('\n'.join(complex_model +
                           [f"LF {ofile2.replace('_params.txt','')}",
                            'Run']))
        
    cmd1 = f"{bt_exe} {intree} {metadata} < {ofile1}"
    if not exists(ofile1.replace('_params.txt', '')+'.Stones.txt'):
        params.append(cmd1)
    cmd2 = f"{bt_exe} {intree} {metadata} < {ofile2}"
    if not exists(ofile2.replace('_params.txt', '')+'.Stones.txt'):
        params.append(cmd2)    
    
with mp.Pool(processes=10) as tp:
    r = list(tqdm(tp.imap(run, params), total=len(params)))


for name, node_str in groups_dict.items():
    ofile1 = join(odir, f'test_{name}_params_{Force_state1}.txt').replace('_params.txt', '')+'.Stones.txt'
    ofile2 = join(odir, f'test_{name}_params_{Force_state2}.txt').replace('_params.txt', '')+'.Stones.txt'
                  
    # larger of the BF, more significant the latter model
    text = summaized_r(complex_f=ofile1,
                       simple_f=ofile2,
                       key='Stone',
                       return_BF_only=True)
    print(name, text)



