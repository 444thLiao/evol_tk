from dating_workflow.toolkit.mcmctree_for import get_posterior_df
import plotly.graph_objects as go
from os.path import join
# p1 = "./dating/scheme1_v2/nucl/clock2_diff_cal/repeat_scheme1_v2_set17_run1"
# p2 = "./dating/scheme1_v2/nucl/clock2_diff_cal/repeat_scheme1_v2_set17_run2_p2"

def compare_two_set(indir1,indir2):
    CI_1 = get_posterior_df(join(indir1, 'mcmc.txt'))
    CI_2 = get_posterior_df(join(indir2, 'mcmc.txt'))

    # remove lnL row
    CI_1 = CI_1.iloc[:-1, :]
    CI_2 = CI_2.iloc[:-1, :]

    dis1 = list(CI_1['Posterior mean time (100 Ma)'])
    dis2 = list(CI_2['Posterior mean time (100 Ma)'])

    fig = go.Figure()
    fig.add_scatter(x=dis1,
                    y=dis2,
                    name='compared',
                    mode='markers')
    fig.add_scatter(x=[min(dis1 + dis2), max(dis1 + dis2)],
                    y=[min(dis1 + dis2), max(dis1 + dis2)],
                    mode='lines',
                    name='y=x')

    fig.layout.width = 1000
    fig.layout.height = 1000
    fig.layout.xaxis.title = "run2 posterior mean time (100Ma)"
    fig.layout.yaxis.title = "run1 posterior mean time (100Ma)"
    
    
    return fig

