import pandas as pd
from ete3 import Tree
from glob import glob
from os.path import *



def parse_mcmc(infile):
    t = pd.read_csv(infile,sep='\t',index_col=0)
    return t.mean()

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
        if row.startswith('t_n'):
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
indir = './dating_for/clock3'
for each_dir in glob(join(indir,'repeat_82g*')):
    if exists(join(each_dir,'FigTree.tre')):
        outfile = glob(join(each_dir,'*.out'))[0]
        mcmc = join(each_dir,'mcmc.txt')
        log = join(each_dir,'run.log')
        t = get_node_name(outfile)
        name = 't_n%s' % t.get_common_ancestor(ns).name
        df = pd.read_csv(mcmc,sep='\t',index_col=0)
        set_name = basename(dirname(outfile)).partition('_')[-1]
        collect_[set_name]= df.loc[:,name]
        df = get_CI(log)
        tmp_df.loc[0,set_name] = '%s (%s) '% (df.loc[name,'Posterior mean time (100 Ma)'],df.loc[name,'CIs'])
    elif exists(join(each_dir,'mcmc.txt')):
        mcmc = join(each_dir,'mcmc.txt')
        name = 't_n217' 
        df = pd.read_csv(mcmc,sep='\t',index_col=0)
        set_name = basename(dirname(mcmc)).partition('_')[-1]
        unfinished_c[set_name]= df.loc[:,name]

tmp_df = tmp_df.T
tmp_df = tmp_df.reindex([_ for _ in tmp_df.index if 'run1' in _])
tmp_df.index = [_.split('_')[1] for _ in tmp_df.index]
tmp_df = tmp_df.reindex(sorted(tmp_df.index,key=lambda x: int(x.split('set')[-1])))
tmp_df.columns = ["divergence time/100Mya (CI)"]
tmp_df.to_excel('./test.xlsx')

tmp = []
for k,v in collect_.items():
    _df = pd.DataFrame(v.values)
    _df.columns = ['time']
    _df.loc[:,'set'] = k.rpartition('_')[0]
    tmp.append(_df)
df = pd.concat(tmp,axis=0)
df.loc[:,'num_set'] = [int(_.split('_set')[-1]) for _ in df.loc[:,'set']]
df.loc[:,'set'] = [_.split('_')[-1].capitalize() for _ in df.loc[:,'set']]
df = df.sort_values('num_set')
import plotly.express as px
fig = px.violin(df,x='set',y='time', box=True,points=False)
fig.layout.yaxis.title.text = 'Divergence time(100Mya)'
fig.layout.xaxis.title.text = 'Sets of calibration information'
fig.layout.yaxis.title.font.size = 30
fig.layout.xaxis.title.font.size = 30
fig.layout.xaxis.tickfont.size = 20
fig.write_html('./82g.html',include_plotlyjs='cdn')