import itertools
from tqdm import tqdm
import numpy as np
from Bio import AlignIO
from collections import defaultdict
def count_same(s1,s2):
    c = 0
    gap_c = 0
    for b1,b2 in zip(s1.seq,s2.seq):
        if b1==b2 and b1 in 'actgACTG':
            c +=1
        elif b1==b2 and b1 == '-':
            gap_c +=1
    return c,gap_c


def get_pairwise_iden(aln):
    iden_pairwise = defaultdict(dict)
    for _1,_2 in tqdm(itertools.combinations(range(len(aln)),2),
                 total=(len(aln)**2-len(aln))/2):
        s1,s2 = aln[_1],aln[_2]
        num_same,gap_c = count_same(s1,s2)
        iden_pairwise[s1.id][s2.id] = float(num_same/(aln.get_alignment_length()-gap_c) )*100
    return iden_pairwise   

aln_16S = AlignIO.read('/home-user/mxie/rpob/16S_all.aln','fasta')

aln_rpoB = AlignIO.read('/home-user/mxie/rpob/all_rpoB-D.aln','fasta')
aln_rpoB.get_alignment_length()

iden_pairwise_rpoB = get_pairwise_iden(aln_rpoB)
iden_pairwise_16S = get_pairwise_iden(aln_16S)


def get_genome2genome(iden_pairwise):
    genome2genome = defaultdict(lambda : defaultdict(list))
    for l1,_d in tqdm(iden_pairwise.items()):
        for l2,iden in _d.items():
            g1 = l1.split('_')[0]   # for iden which is 'genome_locus num'
            g2 = l2.split('_')[0]
            genome2genome[g1][g2].append(iden)
            genome2genome[g2][g1].append(iden)  
    return genome2genome      
# check genomes with multiple target genes
# make sure that the following average process is appropriate

genome2locus = defaultdict(list)
for l1,_d in tqdm(iden_pairwise_16S.items()):
    g1 = l1.split('_')[0]
    genome2locus[g1].append(l1)
multiple_genomes = [g for g,l in genome2locus.items() if len(l)>1]
print(len(multiple_genomes))

multiple_genomes = set(multiple_genomes)
genomes_in_pairwise = defaultdict(list)
for l1,_d in tqdm(iden_pairwise_16S.items()):
    for l2,iden in _d.items():
        g1 = l1.split('_')[0]
        g2 = l2.split('_')[0]
        if g1==g2 and (g1 in multiple_genomes):
            genomes_in_pairwise[g1].append((l1,l2,iden))


genome2genome_16S = get_genome2genome(iden_pairwise_16S)
genome2genome_16S = {k:{_k:np.mean(_v) 
                        for _k,_v in v.items()} 
                     for k,v in genome2genome_16S.items()}

genome2genome_rpoB = get_genome2genome(iden_pairwise_rpoB)
genome2genome_rpoB = {k:{_k:np.mean(_v) for _k,_v in v.items()} for k,v in genome2genome_rpoB.items()}        


shared_genomes = set(genome2genome_rpoB).intersection(set(genome2genome_16S))
print(len(shared_genomes))


## visualizations   (very project specific)
def assign_color(g1,g2):
    cri = set([g1,g2])
    if cri == {'target'}:
        c = '#558B2F'
    elif 'target' in cri and len(cri)>1:
        c = '#2979FF'
    else:
        c = '#DD2C00'    
    return c

color2lineage = {'#0000ff':'outgroup',
                 '#00ff00':'target',
                '#ff0000':'sister',}
genome2lineage = {genome:color2lineage[color] 
                  for genome,color in genome2color.items()}
   
x = []
y = []
t = []
cs = []
for g1,g2 in itertools.combinations(shared_genomes,2):
    x.append(genome2genome_rpoB[g1][g2])
    y.append(genome2genome_16S[g1][g2])
    t.append(f"{g1} to {g2}")
    c = assign_color(genome2lineage[g1],genome2lineage[g2])
    cs.append(c)

import pandas as pd
import plotly.express as px

df = pd.DataFrame()
df.loc[:,'x'] = x
df.loc[:,'y'] = y
df.loc[:,'text'] = t
df.loc[:,'color'] = cs
sub_df = df.loc[(df['x']>40) & (df['y']>88),: ]
dss3_clade = []
for _ in sub_df['text']:
    t = set(_.split(' to '))
    if len(t.intersection(set(clade_genomes)))==1:
        dss3_clade.append('dss3 clade (within group)')
    elif len(t.intersection(set(clade_genomes)))==2:
        dss3_clade.append('dss3 clade (in group)')
    else:
        dss3_clade.append('Not')

# sub_df = df.copy()
sub_df.loc[:,'is_dss'] = ['contain DSS3' if 'GNM000011965' in _ else 'Not' for _ in sub_df['text'] ]
sub_df.loc[:,'dss3_clade'] = dss3_clade
sub_df = sub_df.loc[sub_df['color']=='#558B2F']
fig = px.scatter(sub_df,x='x',y='y',
                 #color='color',
                 color_discrete_map={k:k for k in sub_df['color'].unique()},
                 color='dss3_clade',
                 #marginal_x='box',
                 marginal_y='violin',
                 opacity=0.6
                )

fig.layout.xaxis.range=[70,100.5]
fig.layout.yaxis.range=[92,100.5]
fig.layout.xaxis.autorange=False
fig.layout.xaxis.autorange=False

# fig.write_html('./ipynb/scatter_test.html')
