from dating_workflow.toolkit.mcmctree_for import get_posterior_df
from glob import glob
from os.path import *

import numpy as np
import pandas as pd
import plotly.figure_factory as ff
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from tqdm import tqdm
from collections import defaultdict


prior_set = '~/data/AOB/dating/160g/batch_prior_nucl/*/mcmc.txt'
posterior_set = '~/data/AOB/dating/160g/nucl/clock2_diff_cal/*run1/mcmc.txt'

prior_df = []
prior_list = glob(expanduser(prior_set))
for mcmc in tqdm(prior_list):
    name = basename(dirname(mcmc)).replace('_prior', '')
    _df = pd.read_csv(mcmc, sep='\t', index_col=0)
    _df = _df.reindex(columns=[_ for _ in _df.columns if _.startswith('t_n')])
    prior_df.append((name, _df))
# prior_df = list(sorted(prior_df,key=lambda x: int(x[0].split('set')[-1])))
prior_df = dict(prior_df)

posterior_CIs = {}
posterior_df = []
f = glob(expanduser(posterior_set))
for mcmc in tqdm(f):
    if exists(join(dirname(mcmc), 'FigTree.tre')):
        name = basename(dirname(mcmc)).replace('_run1', '').replace('repeat_', '')
        _df = pd.read_csv(mcmc, sep='\t', index_col=0)
        _df = _df.reindex(columns=[_ for _ in _df.columns if _.startswith('t_n')])
        posterior_df.append((name, _df))
        posterior_CIs[name] = get_posterior_df(join(dirname(mcmc), 'mcmc.txt')))
# posterior_df = list(sorted(posterior_df,key=lambda x: int(x[0].split('set')[-1])))
posterior_df = dict(posterior_df)

shared_set = set(prior_df).intersection(posterior_df)
shared_set = list(sorted(shared_set, key=lambda x: int(x.split('set')[-1])))


group_info = {"set14": ["set14",'set17', 'set19', 'set20', 'set21','set34', 'set23'],
              "set1":['set1', 'set2', 'set4', 'set5', 'set6', 'set8', 'set9', 'set33'],
              "set24":['set24', 'set26', 'set28', 'set29', 'set30', 'set32', 'set36'],
              "sub":['set33', 'set34', 'set35', 'set36', 'set37']
               }

for k,v in group_info.items():
    group_info[k] = [f"83g_{_}" for _ in v]

for _title,shared_set in group_info.items():
    a = pd.read_excel('./dating_for/calibrations_set/calibrations_sets.xlsx', index_col=0)

    fig = make_subplots(rows=len(shared_set), cols=5,
                        shared_xaxes=True,
                        subplot_titles=['Node1',
                                        'Node2',
                                        'Node3',
                                        'Node4',
                                        'Crown group of Anammox'],
                        vertical_spacing=0.003)
    interested_columns = ['t_n84',
                          't_n85',
                          't_n104',
                          "t_n99",
                          't_n139']
    # interested_columns = ["t_n76",
    #                       "t_n77",
    #                       "t_n84",
    #                       "t_n93",
    #                       "t_n123"]


    xmax_dict = defaultdict(int)
    for s in shared_set:
        for n in interested_columns:
            m1 = posterior_df[s].loc[:, n].max()
            m2 = prior_df[s].loc[:, n].max()
            m = np.max([m1, m2])
            if xmax_dict[n] <= m:
                xmax_dict[n] = m

    ys = ['prior ages', 'posterior ages']
    for idx1, s in tqdm(enumerate(shared_set)):
        for idx2, n in enumerate(interested_columns):
            xs = []
            xs.append(list(prior_df[s].loc[:, n]))
            xs.append(list(posterior_df[s].loc[:, n]))
            _fig = ff.create_distplot(xs,
                                      ys,
                                      show_hist=False,
                                      show_rug=False,
                                      bin_size=1,
                                      )
            fig.append_trace(_fig.data[0], row=idx1 + 1, col=idx2 + 1)
            fig.append_trace(_fig.data[1], row=idx1 + 1, col=idx2 + 1)
            # xmax = np.max([_ for _d in _fig.data for _ in _d.x])
            ymax = np.max([_ for _d in _fig.data for _ in _d.y])

            CI = posterior_CIs[s].loc[n, 'CIs']
            CI = '%s - %s' % tuple(map(lambda x: round(float(x), 2),
                                       CI.split(' - ')))
            if idx2 != 4:
                prior_set = a.loc[s.split('_')[-1], :].values[idx2]
            else:
                prior_set = 'nan'
            fig.append_trace(go.Scatter(x=[xmax_dict[n] + 5],
                                        y=[ymax],
                                        text=f"{CI} <Br> ({prior_set})" if str(prior_set) != 'nan' else f"{CI}",
                                        mode='text',
                                        textposition="bottom left"
                                        ),
                             row=idx1 + 1, col=idx2 + 1)
    for idx, v in enumerate(shared_set):
        fig.get_subplot(idx + 1, 1).yaxis.title.text = v.split('_')[-1]

    fig.layout.showlegend = False
    fig.write_html(f'./dating_for/83g/nucl/clock2_diff_cal/{_title}_compare2prior.html', include_plotlyjs='cdn')

fig.layout.height = 2000
fig.write_html('./dating_for/83g/nucl/clock2_diff_cal/compare2prior.html', include_plotlyjs='cdn')

#
# def main(odir,interested_columns,shared_set):
#
#     if not exists(odir):
#         os.makedirs(odir)
#
#     #interested_columns = 't_n84'
#     ys = ['prior ages','posterior ages']
#     for s in tqdm(shared_set):
#         xs = []
#         xs.append(list(prior_df[s].loc[:,interested_columns]))
#         xs.append(list(posterior_df[s].loc[:,interested_columns]))
#
#         fig = ff.create_distplot(xs,
#                                 ys,
#                                 show_hist=False,
#                                 show_rug=False,
#                                 bin_size=1)
#         fig.write_image(join(odir,f'{s}_compare.png'))

if __name__ == "__main__":
    odir = expanduser('~/data/plancto/dating_for/83g/comparison_prior_posterior/clock3_root')
    interested_columns = 't_n84'
    # main(odir)
