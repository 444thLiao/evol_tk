import pandas as pd
from glob import glob
import plotly.express as px
import plotly.graph_objects as go

def get_CI(f):
    f = open(f).read().split('\n')
    head = 'Posterior means (95% Equal-tail CI) (95% HPD CI) HPD-CI-width'
    if head not in f:
        print('no head')
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

def get_r2(text):
    r_squre_text = text['hovertemplate'].split('<br>')[2]
    coef = text['hovertemplate'].split('<br>')[1]
    coef = coef.replace('CI_width','y')
    coef = coef.replace('Posterior mean time (100 Ma)', 'x').split('+')[0].strip()
    prefix, v = r_squre_text.split('=')

    return f"{coef}, {prefix} = {v}"


def draw_r(df):
    df.loc[:,'Posterior mean time (100 Ma)'] = df.loc[:,'Posterior mean time (100 Ma)']/10
    df.loc[:, 'CI_width'] = df.loc[:, 'CI_width'].astype(float) / 10
    fig = px.scatter(df,
                     x='Posterior mean time (100 Ma)',
                     y='CI_width',
                     trendline="ols")

    r_squre_text = [get_r2(_)
                        for _ in [fig.data[1]]]

    fig.update_layout(
        showlegend=False,
        annotations=[
            go.layout.Annotation(
                x=0.5,
                y=max(list(df.loc[:, 'CI_width'])),
                text='<br>'.join(r_squre_text),
                showarrow=False,
                font=dict(size=15)
            )
        ]
    )
    fig.layout.xaxis.title.text = 'Mean posterior divergence time (Gya)'
    fig.layout.yaxis.title.text = '95% HPD CI width (Gyr)'
    r_squre_v = float(r_squre_text[0].split('=')[-1].strip())
    return fig,r_squre_v

from os.path import *
import os
odir = './dating_for/result_infinite_plot'
if not exists(odir):
    os.makedirs(odir)

tmp_df = pd.DataFrame()
a = glob("./dating_for/clock3/*_run1/run.log")
for f1 in a:
    if '187' not in f1:
        name = f1.split('_')[2]
        set_name = f1.split('_')[3]
        f2 = f1.replace('_run1','_run2')
        df1,df2 = get_CI(f1),get_CI(f2)
        time1,time2 = df1.iloc[:,-1],df2.iloc[:,-1]
        fig,r_squre_v = draw_r(df1)
        fig.write_image(join(odir,f'repeat_{name}_{set_name}.png'))
        tmp_df.loc[name,set_name] = r_squre_v

tmp_df = tmp_df.reindex(columns=sorted(tmp_df.columns,key=lambda x:int(x.replace('set',''))))
tmp_df.to_excel('./infinite_site_r2.xlsx')