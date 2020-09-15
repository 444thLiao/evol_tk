from dating_workflow.bin.parse_mcmc import get_CI
import plotly.graph_objects as go
from os.path import join
# p1 = "./dating/scheme1_v2/nucl/clock2_diff_cal/repeat_scheme1_v2_set17_run1"
# p2 = "./dating/scheme1_v2/nucl/clock2_diff_cal/repeat_scheme1_v2_set17_run2_p2"

def compare_two_set(indir1,indir2):
    CI_1 = get_CI(join(indir1, '03_mcmctree.out'))
    CI_2 = get_CI(join(indir2, '03_mcmctree.out'))

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
                    mode='lines')


    return fig