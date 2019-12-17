import pandas as pd
from ete3 import Tree

from global_search.classification_script import diff_marine_non_marine

metadata = '../rawdata/biometadata_952.xlsx'
intree = '../trees/iqtree/over20p_bac120.formatted.newick'
output_txt = './m2nm_{mode}.txt'

metadata_df = pd.read_excel(metadata, index_col=None)
new_df = diff_marine_non_marine(metadata_df)
id2habitat = dict(zip(new_df.iloc[:, 1],
                      new_df.loc[:, 'habitat(auto,diff marine/non-marine)']
                      )
                  )

id2habitats = {k: [v] if not pd.isna(v) else ['unknown'] for k, v in id2habitat.items()}

# outgroup
id2habitat.update({"GCA_000019965.1": 'non-marine',
                   "GCA_000020225.1": 'non-marine',
                   "GCA_000172155.1": 'non-marine',
                   "GCA_001318295.1": 'non-marine',
                   "GCA_001613545.1": 'non-marine',
                   "GCA_900097105.1": 'non-marine',
                   "GCA_001746835.1": 'non-marine'})

# filled multistates.txt
rows = []
tree = Tree(intree, format=3)
for l in tree.get_leaves():
    name = l.name
    habitat = id2habitat.get(name, '-')
    if habitat == 'marine':
        rows.append('\t'.join([name, 'M']))
    elif habitat == 'non-marine':
        rows.append('\t'.join([name, 'N']))
    else:
        rows.append('\t'.join([name, 'NM']))
with open(output_txt.format(mode='multistate'), 'w') as f1:
    f1.write('\n'.join(rows))

# filled multistates.txt
rows = []
tree = Tree(intree, format=3)
for l in tree.get_leaves():
    name = l.name
    habitat = id2habitat.get(name, '-')
    if habitat == 'marine':
        rows.append('\t'.join([name, '1', '0']))
    elif habitat == 'non-marine':
        rows.append('\t'.join([name, '0', '1']))
    else:
        rows.append('\t'.join([name, '-', '-']))
with open(output_txt.format(mode='binary'), 'w') as f1:
    f1.write('\n'.join(rows))
