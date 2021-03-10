"""
it could read in both prior
"""
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
from dating_workflow.figtree2itol import get_node_name
import click

prior_set = '~/data/AOB/dating/160g/batch_prior_nucl/*/mcmc.txt'
posterior_set = '~/data/AOB/dating/160g/nucl/clock2_diff_cal/*run1/mcmc.txt'


def read_multi_mcmc(f_pattern,
                    rename_f= lambda x: basename(dirname(x)).replace('_prior', '').split('_')[-1]):
    df_list = []
    CIs_dict = {}
    mcmc_list = glob(expanduser(f_pattern))
    for mcmc in tqdm(mcmc_list):
        name = rename_f(mcmc)
        _df = pd.read_csv(mcmc, sep='\t', index_col=0)
        _df = _df.reindex(columns=[_ for _ in _df.columns if _.startswith('t_n')])
        df_list.append((name, _df))
    # prior_df = list(sorted(prior_df,key=lambda x: int(x[0].split('set')[-1])))
        CIs_dict[name] = get_posterior_df(mcmc)
    df_dict = dict(df_list)
    return df_dict,CIs_dict
prior_df,prior_CIs = read_multi_mcmc(prior_set)
posterior_df,posterior_CIs = read_multi_mcmc(posterior_set,
                           rename_f=lambda x: basename(dirname(x)).replace('_run1', '').replace('repeat_', '').split('_')[-1])

shared_set = set(prior_df).intersection(posterior_df)
shared_set = list(sorted(shared_set, key=lambda x: int(x.split('set')[-1])))
shared_set = [_.split('_')[-1] for _ in shared_set]


group_info = {"set1":  ["set1"] + [f"set{n}" for n in range(5,16)],
              "set2": ["set2"] + [f"set{n}" for n in range(16,26)],
               }

# for k,v in group_info.items():
#     group_info[k] = [f"160g_{_}" for _ in v]
_interested_columns = ['GCA_000011385.1|GCA_004357655.1',
                       'GCA_000011385.1|GCA_000013205.1',
                       'GCA_000196515.1|GCA_001548455.1',
                       "GCA_000317575.1|GCA_000317025.1",
                       'GCA_003695135.1|GCA_900172325.1',
                       "GCA_900172325.1|GCA_002215215.1"]
interested_columns = []
t = get_node_name(f'./dating/160g/nucl/clock2_diff_cal/repeat_160g_set10_run1/03_mcmctree.out')
for _ in _interested_columns:
    name = 't_n%s' % t.get_common_ancestor(_.split('|')).name
    interested_columns.append(name)

for _title,shared_set in group_info.items():
    a = pd.read_excel('./dating/calibrations_set/calibrations_sets.xlsx', index_col=0)
    fig = make_subplots(rows=len(shared_set), cols=len(interested_columns),
                        shared_xaxes=True,
                        subplot_titles=['Root',
                                        'Node1',
                                        'Node2',
                                        'Node3',
                                        'Node4',
                                        'Node5',],
                        vertical_spacing=0.003)

    xmax_dict = defaultdict(int)
    for set_name in shared_set:
        for n in interested_columns:
            m1 = posterior_df[set_name].loc[:, n].max()
            m2 = prior_df[set_name].loc[:, n].max()
            m = np.max([m1, m2])
            if xmax_dict[n] <= m:
                xmax_dict[n] = m

    ys = ['prior ages', 'posterior ages']
    for idx1, set_name in tqdm(enumerate(shared_set)):
        for idx2, n in enumerate(interested_columns):
            prior_set = a.loc[set_name.split('_')[-1], :].values[idx2]
            # don't draw distribution of the node without constraints
            if str(prior_set) == 'nan':
                continue

            xs = []
            xs.append(list(prior_df[set_name].loc[:, n]))
            xs.append(list(posterior_df[set_name].loc[:, n]))
            _fig = ff.create_distplot(xs,
                                      ys,
                                      show_hist=False,
                                      show_rug=False,
                                      bin_size=1,
                                      )

            ymax = np.max([_ for _d in _fig.data for _ in _d.y])
            CI = posterior_CIs[set_name].loc[n, 'CIs']
            CI = '%s - %s' % tuple(map(lambda x: round(float(x), 2),
                                       CI.split(' - ')))

            fig.append_trace(_fig.data[0], row=idx1 + 1, col=idx2 + 1)
            fig.append_trace(_fig.data[1], row=idx1 + 1, col=idx2 + 1)
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
    fig.layout.height = 2000
    fig.write_html(f'./dating/160g/nucl/clock2_diff_cal/{_title}_compare2prior.html', include_plotlyjs='cdn')


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
