"""
Raw script
Project specific

"""

from collections import defaultdict
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq
import re, itertools
from tqdm import tqdm
from Bio import AlignIO
import pandas as pd
from skbio import DNA, TabularMSA

records = list(SeqIO.parse("all_rpoB-D.aln", "fasta"))
msa = TabularMSA([DNA(str(_.seq).upper()) for _ in records])
msa.index = [_.id for _ in records]

df = pd.read_excel("./table01,808_genome_summry,20191211,XYF.xlsx")
df = df.set_index("GNMUID")
all_faa = list(SeqIO.parse("./all_rpoB-D.fas", "fasta"))
seq2target = {_.id: df.loc[_.id.split("_")[0], "tax"] for _ in all_faa}
# text = to_color_strip(seq2target,
#                       info2color={'sister group':"#ff0000","target":"#00ff00","outgroup":"#0000ff"})
# with open('./colorstrip_gene.txt','w') as f1:
#     f1.write(text)
all_group = set(seq2target.values())
g2seqs = {}
for g in all_group:
    group = [_ for _, v in seq2target.items() if v == g]
    g2seqs[g] = group

frame_list = list(zip(range(msa.shape[1]), range(10, msa.shape[1])))
x = []
y_d = defaultdict(list)
for s, e in tqdm(frame_list):
    x.append((s + e) / 2)
    # try:
    for name, group in g2seqs.items():
        con = msa.loc[group, s:e].conservation(gap_mode="ignore").mean()
        y_d[name].append(con)
    for g1, g2 in list(itertools.combinations(g2seqs, 2)):
        con = (
            msa.loc[g2seqs[g1] + g2seqs[g2], s:e].conservation(gap_mode="ignore").mean()
        )
        y_d[f"{g1}+{g2}"].append(con)


import plotly.graph_objs as go

fig = go.Figure()
for name in y_d:
    fig.add_scatter(x=x, y=y_d[name], name=name)
fig.show()
