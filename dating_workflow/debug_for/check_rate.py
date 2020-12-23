"""
For debugging the inner mechanism of MCMCTree
Here is a script for investigating the different mutation rates across branches.


"""
import pandas as pd 
import plotly
import plotly.graph_objects as go
import plotly.figure_factory as ff

import numpy as np
from scipy import stats

f = './03_mcmctree.log'
def get_rates_df(f):
    f = open(f).read().split('\n')
    head = 'Posterior means (95% Equal-tail CI) (95% HPD CI) HPD-CI-width'
    idx = f.index(head)
    remained_txt = f[idx + 1:]

    def format_v(x):
        x = x.strip('(')
        x = x.strip(')')
        x = x.strip(',')
        return float(x)

    idx = []
    CIs = []
    mean_collect = []
    CI_collect = []
    for row in remained_txt:
        if row.startswith('t_n') or row.startswith('lnL') or row.startswith('r_n') or  row.startswith('mu') or  row.startswith('sigma2'):
            vals = row.split(' ')
            vals = [_ for _ in vals if _ and _ not in '(),']
            posterior_mean, equal_tailed_95p_down, equal_tailed_95p_up = map(format_v, vals[1:4])
            CI_collect.append((equal_tailed_95p_up - equal_tailed_95p_down))
            mean_collect.append(posterior_mean)
            idx.append(vals[0])
            CIs.append('%s - %s' % (equal_tailed_95p_down, equal_tailed_95p_up))
    df = pd.DataFrame()
    df.loc[:, 'CI_width'] = CI_collect
    df.loc[:, 'CIs'] = CIs
    df.loc[:, 'Posterior mean time (100 Ma)'] = mean_collect
    df.index = idx

    time_df = df.reindex([_ for _ in df.index if _.startswith('t_n')])
    rates_df = df.reindex([_ for _ in df.index if _.startswith('r_n') ]) # no root
    return rates_df 


rates_df = get_rates_df('../prior/nodata_mcmctree.log')
# y = list(-np.log(rates_df.iloc[:,2]))
y = list(rates_df.iloc[:,2])
fig = ff.create_distplot([y], ['rate'], 
                         show_hist=False, )
fig.show()


y = -np.log(rates_df.iloc[:,2])
shape, loc, scale = stats.lognorm.fit(y, floc=0)
mu = np.log(scale)
sigma = shape

print(sigma,mu)
