"""
Python version of ClonalFrameML/src/cfml_results.R

is it needed?

"""

from Bio import SeqIO
import pandas as pd
import numpy as np
from ete3 import Tree

xref = np.array(open("MC15.position_cross_reference.txt").read().strip().split(',')).astype(int) -1 # make the 1-coordinates into 0-coordinates
state_df = pd.read_csv("./MC15.importation_status.txt",sep='\t')
records = list(SeqIO.parse('./MC15.ML_sequence.fasta','fasta'))
seq_list = [list(str(_.seq)) for _ in records]
seq_array = np.array(seq_list)
seq_df = pd.DataFrame(seq_array,index=[_.id for _ in records])

tre = Tree('MC15.labelled_tree.newick',1)
n2node = {n.name:n for n in tre.traverse()}
node2parent = [n2node[n].up.name if n2node[n].up else seq_df.index[0] for n in seq_df.index]
# note that the seq_df.index[0] is randomly given and it is meanless Since it will be replaced as 0 with the following lines
# whether mutation
wh_mut = seq_df.apply(lambda x:x.values!=x[node2parent].values,axis=0).astype(int)
wh_mut.loc[seq_df.index[-1],:] = 0 # mark the last row as 0

# number of mutation (might be 0)
n_mut = wh_mut.loc[[_ for _ in wh_mut.index if _.startswith('NODE_')],:].sum()  # only Sum the mutation at the internal nodes
is_homoplasy = n_mut>1

# spectrum of mutation that are masked into 0 if it is not a mutation (from whether mutations)
pos_with_diffbp = xref[xref>-1]  
a = wh_mut.iloc[:,pos_with_diffbp]  # each node versus each real positions (might be duplicated)
b = n_mut[pos_with_diffbp]  # each node versus each real positions (might be duplicated)
spectrum_mut = a*b

# 
from visualization.tanglegram import tree_vis
import plotly.graph_objects as go
tv = tree_vis(tre)
data = tv.get_plotly_data()
fig = go.Figure()
fig.add_traces(data)

