
from os.path import *
import pandas as pd
from dating_workflow.toolkit.mcmctree_for import get_node_name_from_log,read_mcmc
import plotly.express as px

def get_fig(_df):
    _df = _df.loc[_df["group name"].isin(list(g2color)), :]
    _fig = px.violin(
        _df,
        y="cal",
        x="time",
        color="group name",
        color_discrete_map=g2color,
        points=False,
        orientation="h"
        # color='calibration sets',
        # showlegend=False
    )
    num_y = len(_df["cal"].unique())
    _fig.update_traces(
        side="positive",
        width=1.8,
    )
    return _fig

indir = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/dating/sys_testing_latest/MITO24/B12E5.euk_cyano_IR_nucl'
tre = get_node_name_from_log(f'{indir}/prior/nodata_mcmctree.log')
dfs = []
for cal in [f"B{i}E{s}" for i in range(1,15) for s in range(1,7) if i not in [9,10]
           ]:
    mcmc = f"{dirname(indir)}/{cal}.euk_cyano_IR_nucl/mcmc_for/mcmc.txt"
    if not exists(mcmc):continue
    df = read_mcmc(mcmc)
    df = df.sample(5000)
    
    for lca,name in [#('GCA_001828545.1,GCA_005524015.1','Anammox'),
                     ('GCA_013697045.1,GCA_001644685.1','Gamma-AOB'),
                     ('GCA_001772005.1,GCA_013521015.1','Beta-AOB'),
                     ('GCA_017879665.1,GCA_013140535.1','Comammox'),
                    ]:
        node = tre.get_common_ancestor(lca.split(','))
        #print(name,node.name,len(node.get_leaf_names()))
        times = df[[f"t_n{node.name}"]]
        times.columns = ['time']
        times.loc[:,'group name'] = name
        times.loc[:,'cal'] = cal
        dfs.append(times)
_df = pd.concat(dfs,axis=0)
g2color = {
    "Gamma-AOB": "#78fce0",
    "Beta-AOB": "#956bb4",
    "Comammox": "#edc21a",
    #"Anammox": "#ff8000",
}

full_fig = get_fig(_df,)
full_fig.add_scatter(x=[27,27],y=['B1E1','B14E6'],mode='lines',opacity=0.3,line=dict(color='#000000'))   
full_fig.layout.template = "simple_white"
#full_fig.layout.xaxis.range = [60, 0]
full_fig.layout.xaxis.autorange = 'reversed'
full_fig.layout.width = 700
full_fig.layout.height = 650
full_fig.show()