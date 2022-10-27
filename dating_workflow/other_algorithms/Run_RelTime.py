"""
Run RelTime implemented in MEGA X.

Not executable scripts. Just as reference
Use MEGA-GC (command-line interface)

"""


import os
from ete3 import Tree
from tqdm import tqdm
import re

ref_tree = Tree("/mnt/home-backup/thliao/wol_db/data/trees/astral/branch_length/cons/astral.cons.nid.nwk",1)
n = [_ for _ in ref_tree.traverse() if _.name == 'N7'][0]  # LCA of bac
assert len(n.get_leaf_names()) == 8452

outgroups = n.up.children[0].get_leaf_names()
with open('./tmp.txt','w') as f1:
    for k in outgroups:
        f1.write(f"{k}=outgroup\n")

calibration = f"""!MRCA='Root' TaxonA='G000212395' TaxonB='G001791795' maxtime=45.00000000 calibrationName='Root-split';
!MRCA='Total cyano' TaxonA='G000011385' TaxonB='G001786505' Distribution=uniform mintime=23.00000000 maxtime=38.00000000 calibrationName='Total cyano-split';
!MRCA='Total Nostoc' TaxonA='G000952155' TaxonB='G000734895' Distribution=uniform mintime=16.00000000 maxtime=38.00000000 calibrationName='Nostoc-split';
!MRCA='Total Pleuro' TaxonA='G000317575' TaxonB='G000317655' Distribution=uniform mintime=17.00000000 maxtime=38.00000000 calibrationName='Pleuro-split';
"""
with open('./wol_cal.txt','w') as f1:
    f1.write(calibration)

# for g in re.findall('G[0-9]+',calibration):
#     if g in outgroups:
#         print(g)

# conversion (since the tree required by MEGA can not contain the internal node name, we need to convert the generated time tree.)
tre = Tree(f'/mnt/ivy/thliao/project/ML_oxygen/testing_sets/WoL/time_wol_named.nwk',3)
nt = Tree(f'/mnt/ivy/thliao/project/ML_oxygen/testing_sets/WoL/root45_exactTimes.nwk',quoted_node_names=1)
for n in tqdm(nt.traverse()):
    if not n.is_leaf():
        n2 = tre.get_common_ancestor(n.get_leaf_names())
        assert set(n2.get_leaf_names()) == set(n.get_leaf_names())
        n.name = n2.name
nt.write(outfile=f'/mnt/ivy/thliao/project/ML_oxygen/testing_sets/WoL/root45_dated.nwk',format=3)


reltime_nt = Tree(f'/mnt/ivy/thliao/project/ML_oxygen/testing_sets/WoL/root45_relTimes.nwk',quoted_node_names=1)
for n in tqdm(reltime_nt.traverse()):
    if not n.is_leaf():
        n2 = tre.get_common_ancestor(n.get_leaf_names())
        assert set(n2.get_leaf_names()) == set(n.get_leaf_names())
        n.name = n2.name
reltime_nt.write(outfile=f'/mnt/ivy/thliao/project/ML_oxygen/testing_sets/WoL/root45_reltime_renamed.nwk',format=3)


