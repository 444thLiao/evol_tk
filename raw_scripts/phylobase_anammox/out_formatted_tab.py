import pandas as pd
from ete3 import Tree
from glob import glob
import os




def get_dates(chain_name):
    label_tree = Tree(f"{chain_name}_sample.labels",format=1)
    dates = pd.read_csv(f"{chain_name}_sample.dates",sep='\t',index_col=0)
    # rate_tree = Tree(f"{chain_name}_sample.ratetree",format=1)


    name2group = {"anammox bacteria": "GCA_001828545.1|GCA_004282745.1",
                  "root":"GCA_000011385.1|GCA_003576915.1",
                  "cyanobacteria":"GCA_000011385.1|GCA_000013205.1",
                  "Nostocales":"GCA_000196515.1|GCA_001548455.1",
                  "pleurocapsales":"GCA_000317575.1|GCA_000317025.1",}

    c = []
    for gname, group in name2group.items():
        group = group.split('|')
        raw_name = '%s' % label_tree.get_common_ancestor(group).name
        sub_dates = dates.loc[[raw_name],:]
        sub_dates.index = [gname]
        c.append(sub_dates)

    df = pd.concat(c,axis=0)
    return df
chain_name = "run1/test"
df1 = get_dates(chain_name)
chain_name = "run2/test"
df2 = get_dates(chain_name)
chain_name = "prior_only/test"
prior_df = get_dates(chain_name)

