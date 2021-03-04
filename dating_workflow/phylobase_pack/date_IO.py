

from dating_workflow.toolkit.mcmctree_for import get_node_name_from_log
from ete3 import Tree
import pandas as pd
from os.path import *


# example
chain_name = "run1/test"
date_file = chain_name + '_sample.dates'
label_file = chain_name + "_sample.labels"

t = Tree(label_file,format=1)
in2leaves = {n: n.get_leaf_names() for n in t.traverse() if not n.is_leaf()}

time_estimate = pd.read_csv(date_file,sep='\t',index_col=0)

time = time_estimate["(meandate"]

# inf95, lower bound of 95% Credibility Intervals (CIs)
# sup95, upper bound of 95% Credibility Intervals (CIs)

def get_time(chain_name):
    date_file = chain_name + '_sample.dates'
    label_file = chain_name + "_sample.labels"

    t = Tree(label_file,format=1)
    in2leaves = {n: n.get_leaf_names() for n in t.traverse() if not n.is_leaf()}
    time_estimate = pd.read_csv(date_file,sep='\t',index_col=0)

    time = time_estimate["(meandate"]
    return time,in2leaves
t1,i2l_1 = get_time("run1/test")
t2,i2l_2 = get_time("run2/test")
t_prior,i2l_prior = get_time("prior_only/test")

dates_df = pd.DataFrame()

from dating_workflow.bin.parse_mcmc import get_posterior_df
mcmc = '/mnt/home-backup/thliao/plancto/dating_for/83g/nucl/clock2_diff_cal/repeat_83g_set14_run1/mcmc.txt'
df = get_posterior_df(mcmc)
t = get_node_name_from_log(mcmc.replace('mcmc.txt','run.log'))
n1_2_n2 = {}
for n,l in i2l_1.items():
    _n = t.get_common_ancestor(l)
    n1_2_n2['t_n'+str(_n.name)] = n.name

order_date_phylobayes = t1.reindex([n1_2_n2[_] for _ in df.index if _ in n1_2_n2])
sub_df = df.reindex([_ for _ in df.index if _ in n1_2_n2])

import plotly.express as px
def draw_conv_plot(time1,time2):
    df = pd.DataFrame()
    df.loc[:, 'time1'] = time1/10
    df.loc[:, 'time2'] = time2/10
    fig = px.scatter(df,
                     x='time1',
                     y='time2',
                     #trendline="ols"
                     )
    fig.add_scatter(x=df['time1'],
                    y=df['time1'],
                    mode='lines')
    
    # r_squre_text = [get_r2(_)
    #                     for _ in [fig.data[1]]]

    fig.update_layout(
        showlegend=False,
        # annotations=[
        #     go.layout.Annotation(
        #         x=0.5,
        #         y=max(list(df.loc[:, 'time1'])),
        #         text='<br>'.join(r_squre_text),
        #         showarrow=False,
        #         font=dict(size=15)
        #     )
        # ]
    )
    return fig

fig = draw_conv_plot(time1=sub_df['Posterior mean time (100 Ma)'].values,
             time2 = order_date_phylobayes.values)
