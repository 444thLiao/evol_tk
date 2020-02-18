import pandas as pd
from ete3 import Tree
from glob import glob
from os.path import *
from tqdm import tqdm
import numpy as np
def cal_ESS(df):

    N = int(df.shape[0] * 0.1)
    df = df.iloc[N:, :]
    N = df.shape[0]

    def f(t1):
        def cal(x,k):
            return np.corrcoef(x[k:],x[:-k])[0][1]
        V = np.array([cal(t1.values,k) for k in range(1,N-1)])
        ess = N / (1+sum(V))
        return ess
    ess_array = df.apply(f,axis=0)

def get_CI(f):
    f = open(f).read().split('\n')
    head = 'Posterior means (95% Equal-tail CI) (95% HPD CI) HPD-CI-width'
    if head not in f:
        return None
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
        if row.startswith('t_n') or row.startswith('lnL'):
            vals = row.split(' ')
            vals = [_ for _ in vals if _ and _ not in '(),']
            posterior_mean, equal_tailed_95p_down, equal_tailed_95p_up = map(format_v, vals[1:4])
            CI_collect.append((equal_tailed_95p_up - equal_tailed_95p_down))
            mean_collect.append(posterior_mean)
            idx.append(vals[0])
            CIs.append('%s - %s' % (equal_tailed_95p_down,equal_tailed_95p_up))
    df = pd.DataFrame()
    df.loc[:, 'CI_width'] = CI_collect
    df.loc[:,'CIs'] = CIs
    df.loc[:, 'Posterior mean time (100 Ma)'] = mean_collect
    df.index = idx
    return df

def get_node_name(f):
    matched_row = ''
    match = False
    for row in open(f):
        row = row.strip('\n').strip()
        if match and not row:
            break
        if row.startswith('Species tree for FigTree.  Branch lengths = posterior mean times; 95% CIs = labels'):
            match = True
            continue
        if match:
            matched_row = row
    t = Tree(matched_row.replace(' ', ''), format=8)
    for l in t.get_leaves():
        l.name = l.name.partition('_')[-1]
    return t

tmp_df = pd.DataFrame()
collect_ = {}
unfinished_c = {}
ns = ['GCA_001828545.1','GCA_004282745.1']
indir = './dating_for/83g/clock2_diff_cal/'
interate_dir = list(glob(join(indir,'*'))) #+ ['./dating_for/83g/clock2_rgene/',
                                           #   './dating_for/83g/clock2_sigma2/']
# name = 't_n139'
name = ''
for each_dir in tqdm(interate_dir):
    outfile = glob(join(each_dir, '*.out'))[0]
    mcmc = join(each_dir, 'mcmc.txt')
    log = join(each_dir, 'run.log')
    if exists(join(each_dir,'FigTree.tre')):
        if not name:
            t = get_node_name(outfile)
            name = 't_n%s' % t.get_common_ancestor(ns).name
        df = pd.read_csv(mcmc,sep='\t',index_col=0)
        set_name = basename(dirname(outfile)).partition('_')[-1]
        if not set_name:
            set_name = basename(dirname(outfile))
        collect_[set_name]= df.loc[:,name]
        df = get_CI(log)
        tmp_df.loc[set_name,'Anammox group'] = '%s (%s) '% (df.loc[name,'Posterior mean time (100 Ma)'],
                                                            df.loc[name,'CIs'])
        tmp_df.loc[set_name,'ROOT'] = '%s (%s) '% (df.ix[0,'Posterior mean time (100 Ma)'],
                                                   df.ix[0,'CIs'])
        tmp_df.loc[set_name,'lnL'] = '%s (%s) '% (df.loc["lnL",'Posterior mean time (100 Ma)'],
                                                   df.loc["lnL",'CIs'])
    # if exists(mcmc):
    #     df = pd.read_csv(mcmc,sep='\t',index_col=0)
    #     set_name = basename(dirname(mcmc)).partition('_')[-1]
        # unfinished_c[set_name]= df.loc[:,name]
tmp_df.loc[:,'num_set'] = [int(_.split('_')[1].replace('set',''))
                           if 'run' in _ else 0
                        for _ in tmp_df.index]
tmp_df = tmp_df.sort_values('num_set')
tmp_df.index = [_.replace('83g','83g_clock2')
                for _ in tmp_df.index]
#tmp_df.columns = ["divergence time/100Mya (CI)"]
tmp_df.to_excel('./dating_for/83g/83g_clock2_diff_cal.xlsx',index_label='Calibration set')







tmp = []
for k,v in collect_.items():
    _df = pd.DataFrame(v.values)
    _df.columns = ['time']
    _df.loc[:,'name'] = k
    if 'run1' in k or 'run' not in k:
        tmp.append(_df)
df = pd.concat(tmp,axis=0)
df.loc[:,'num_set'] = [int(_.split('_set')[-1].split('_')[0])
                       if 'run' in _ else 0
                       for _ in df.loc[:,'name']]
df.loc[:,'set'] = [_.replace('83g_','').replace('_run1','').replace('clock2_','') for _ in df.loc[:,'name']]
df = df.sort_values('num_set')

import plotly.express as px
fig = px.violin(df,x='set',y='time', box=True,points=False)
fig.layout.yaxis.title.text = 'Divergence time(100Mya)'
fig.layout.xaxis.title.text = 'Sets of calibration information'
fig.layout.yaxis.title.font.size = 30
fig.layout.xaxis.title.font.size = 30
fig.layout.xaxis.tickfont.size = 20
fig.write_html('./dating_for/83g/83g_clock2_diff_cal.html',include_plotlyjs='cdn')