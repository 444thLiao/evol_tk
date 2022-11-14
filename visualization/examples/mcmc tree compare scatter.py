
from os.path import *
import pandas as pd
from dating_workflow.toolkit.mcmctree_for import get_node_name_from_log,read_mcmc
import plotly.express as px

## need to be adjusted
p1 = ''
p2 = ''
######
tre1 = get_node_name_from_log(f'{p1}/mcmc_for/03_mcmctree.log')
t1 = read_mcmc(f'{p1}/mcmc_for/mcmc.txt')

tre2 = get_node_name_from_log(f'{p2}/mcmc_for/03_mcmctree.log')
t2 = read_mcmc(f'{p2}/mcmc_for/mcmc_for/mcmc.txt')
print(len(tre1.get_leaf_names()),len(tre2.get_leaf_names()))

n2n = {}
for n in tre1.traverse():
    if n.is_leaf():continue
    l,r = n.children
    l = l.get_leaf_names()[0];r=r.get_leaf_names()[0]
    try:
        n2 = tre2.get_common_ancestor([l,r])
        if set(n.get_leaf_names()) == set(n2.get_leaf_names()):
            n2n[n.name] = n2.name
        #print(n.name,n2.name)
    except:
        continue

fig = go.Figure()
xs = [];ys = []
names = [];colors = []
for n,n2 in n2n.items():
    xs.append(t1[f"t_n{n}"].mean())
    ys.append(t2[f"t_n{n2}"].mean())
    names.append(n)
    if n in [str(_) for _ in range(120,129)]:
        # some nodes needed to be highlighted
        colors.append('#f57c00')
    if n in [str(_) for _ in range(135,147)]:
        colors.append('#b2e0df')
    else:
        colors.append('#636efa')
fig.add_scatter(x=xs,y=ys,text=names,mode='markers',marker=dict(color=colors))
fig.add_scatter(x=[0,xs[0]],y=[0,xs[0]],mode='lines')
fig.layout.width=500;fig.layout.height=500
fig.show()