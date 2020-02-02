import os
import warnings
from glob import glob
from os.path import *

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

warnings.filterwarnings('ignore')

infile = './03_mcmctree.out'
ofile = './evaluate.png'


def get_vals(infile):
    f = open(infile).read().split('\n')
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

    mean_collect = []
    CI_collect = []
    for row in remained_txt:
        if row.startswith('t_n'):
            vals = row.split(' ')
            vals = [_ for _ in vals if _ and _ not in '(),']
            posterior_mean, equal_tailed_95p_down, equal_tailed_95p_up = map(format_v, vals[1:4])
            CI_collect.append((equal_tailed_95p_up - equal_tailed_95p_down))
            mean_collect.append(posterior_mean)
    df = pd.DataFrame()
    df.loc[:, 'CI_width'] = CI_collect
    df.loc[:, 'Posterior mean time (100 Ma)'] = mean_collect

    return df


def get_r2(text):
    r_squre_text = text['hovertemplate'].split('<br>')[2]
    prefix, v = r_squre_text.split('=')
    mode_text = text['hovertemplate'].split('<br>')[4]
    mode = mode_text.split('=')[-1]

    return f"{prefix} of {mode}={v}"


def draw_picture(infile, infile2=None, ofile=None, return_fig=False):
    df1 = get_vals(infile)
    df1.loc[:, 'mode'] = 'posterior'
    if infile2 is not None:
        df2 = get_vals(infile2)
        df2.loc[:, 'mode'] = 'prior'
    else:
        df2 = pd.DataFrame()

    df = pd.concat([df1, df2], axis=0)

    fig = px.scatter(df,
                     x='Posterior mean time (100 Ma)',
                     y="CI_width",
                     color="mode",
                     trendline="ols")
    if infile2 is not None:
        r_squre_text = [get_r2(_)
                        for _ in [fig.data[1],
                                  fig.data[3]]]
    else:
        r_squre_text = [get_r2(_)
                        for _ in [fig.data[1]]]
    fig.update_layout(
        showlegend=False,
        annotations=[
            go.layout.Annotation(
                x=0,
                y=max(list(df.loc[:, 'CI_width'])),
                text='<br>'.join(r_squre_text),
                showarrow=False,
                font=dict(size=20)
            )
        ]
    )

    r_squre = float([_.split('=')[-1] for _ in r_squre_text if 'posterior' in _][0])
    if return_fig:
        return fig
    else:
        fig.write_image(ofile)
        return r_squre


collect_r2 = {}
for infile in glob('./locus_each/*/mcmc_for/03_mcmctree.out'):
    infile2 = infile.replace("mcmc_for/03_mcmctree.out", "prior/nodata_mcmctree.out")
    gene_name = infile.split('/')[2].split('_')[-1]
    ofile = f'./pack_result/test_r2/one_gene{gene_name}.png'
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    r2 = draw_picture(infile, infile2, ofile, return_fig=False)
    collect_r2[gene_name] = r2

fs = glob('./design_scheme/*/mcmc_for/03_mcmctree.out')
for infile in fs:
    infile2 = infile.replace("mcmc_for/03_mcmctree.out", "prior/nodata_mcmctree.out")
    gene_name = infile.split('/')[2]
    ofile = f'./pack_result/test_r2/{gene_name}.png'
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    r2 = draw_picture(infile, infile2, ofile, return_fig=False)
    collect_r2[gene_name] = r2

infile = '/home-user/thliao/data/nitrification_for/dating_for/rough_dating/01_mcmctree_out/out'
infile2 = '/home-user/thliao/data/nitrification_for/dating_for/rough_dating/nodata_out/no_data.out'
ofile = './pack_result/test_r2/rough_dating.png'
r2 = draw_picture(infile, infile2, ofile, return_fig=False)
collect_r2['rough dating'] = r2
