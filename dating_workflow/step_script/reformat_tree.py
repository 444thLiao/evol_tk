from ete3 import Tree

tree_file = '243g.formatted.newick'
calibration_txt = '../../calibration.txt'
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

with open('./243g_120gene.calibrations.newick','w') as f1:
    f1.write(f'{len(t.get_leaves())}\n'+final_tree)