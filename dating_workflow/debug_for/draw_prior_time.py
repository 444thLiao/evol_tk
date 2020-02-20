import os
from glob import glob
from os.path import *

import pandas as pd
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
from tqdm import tqdm

prior_df = []
prior_list = glob(expanduser('~/data/plancto/dating_for/83g/batch_prior/*/mcmc.txt'))
for mcmc in tqdm(prior_list):
    name = basename(dirname(mcmc)).replace('_prior','')
    _df = pd.read_csv(mcmc,sep='\t',index_col = 0)
    _df = _df.reindex(columns = [_ for _ in _df.columns if _.startswith('t_n')])
    prior_df.append((name,_df))
prior_df = list(sorted(prior_df,key=lambda x: int(x[0].split('set')[-1])))
prior_df = dict(prior_df)


posterior_df = []
f = glob(expanduser('~/data/plancto/dating_for/83g/clock2_diff_cal/*run1/mcmc.txt'))
for mcmc in tqdm(f):
    if exists(join(dirname(mcmc),'FigTree.tre')):
        name = basename(dirname(mcmc)).replace('_run1','').replace('repeat_','')
        _df = pd.read_csv(mcmc,sep='\t',index_col = 0)
        _df = _df.reindex(columns = [_ for _ in _df.columns if _.startswith('t_n')])
        posterior_df.append((name,_df))
posterior_df = list(sorted(posterior_df,key=lambda x: int(x[0].split('set')[-1])))
posterior_df = dict(posterior_df)


shared_set = set(prior_df).intersection(posterior_df)
shared_set = list(sorted(shared_set,key=lambda x: int(x.split('set')[-1])))


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
ys = ['prior ages','posterior ages']
for idx1,s in tqdm(enumerate(shared_set)):

    for idx2,n in enumerate(interested_columns):
        xs = []
        xs.append(list(prior_df[s].loc[:, n]))
        xs.append(list(posterior_df[s].loc[:, n]))
        _fig = ff.create_distplot(xs,
                                 ys,
                                 show_hist=False,
                                 show_rug=False,
                                 bin_size=1,
                                  )

        fig.append_trace(_fig.data[0],row=idx1+1,col=idx2+1)
        fig.append_trace(_fig.data[1], row=idx1 + 1, col=idx2 + 1)

for idx,v in enumerate(shared_set):
    fig.get_subplot(idx+1,1).yaxis.title.text = v.split('_')[-1]


fig.layout.showlegend = False
fig.layout.height = 2000
fig.write_html('./test.html',include_plotlyjs='cdn')
def main(odir,interested_columns,shared_set):

    if not exists(odir):
        os.makedirs(odir)

    #interested_columns = 't_n84'
    ys = ['prior ages','posterior ages']
    for s in tqdm(shared_set):
        xs = []
        xs.append(list(prior_df[s].loc[:,interested_columns]))
        xs.append(list(posterior_df[s].loc[:,interested_columns]))
        
        fig = ff.create_distplot(xs,
                                ys, 
                                show_hist=False,
                                show_rug=False,
                                bin_size=1)
        fig.write_image(join(odir,f'{s}_compare.png'))

if __name__ == "__main__":
    odir = expanduser('~/data/plancto/dating_for/83g/comparison_prior_posterior/clock3_root')
    interested_columns = 't_n84'
    main(odir)