import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from glob import glob
from tqdm import tqdm
from os.path import *
import os
import warnings
warnings.filterwarnings('ignore')

infile = './03_mcmctree.out'
ofile = './evaluate.png'

def draw_picture(infile,ofile,return_fig=False):
    f = open(infile).read().split('\n')
    head = 'Posterior mean (95% Equal-tail CI) (95% HPD CI) HPD-CI-width'
    if head not in f:
        return None
    idx = f.index(head)
    remained_txt = f[idx+1:]

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
            posterior_mean,equal_tailed_95p_down,equal_tailed_95p_up = map(format_v,vals[1:4])
            CI_collect.append((equal_tailed_95p_up-equal_tailed_95p_down))
            mean_collect.append(posterior_mean)
    
    df = pd.DataFrame()

    df.loc[:,'CI_width'] = CI_collect
    df.loc[:,'Posterior mean time (100 Ma)'] = mean_collect


    fig = px.scatter(df, 
                    x='Posterior mean time (100 Ma)',
                    y="CI_width",
                    trendline="ols")
    r_squre_text = fig.data[1]['hovertemplate'].split('<br>')[2]
    fig.update_layout(
        showlegend=False,
        annotations=[
            go.layout.Annotation(
                x=0,
                y=max(CI_collect),
                text=r_squre_text,
                showarrow=False,
                font=dict(size=20)
            )
        ]
    )
    r_squre = float(r_squre_text.split('=')[-1])
    if return_fig:
        return fig
    else:
        fig.write_image(ofile)
        return r_squre



collect_r2 = {}
for infile in glob('./locus_each/*/mcmc_for/03_mcmctree.out'):
    gene_name = infile.split('/')[2].split('_')[-1]
    ofile = f'./pack_result/test_r2/one_gene{gene_name}.png'
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    r2 = draw_picture(infile,ofile,return_fig=False)
    collect_r2[gene_name] = r2

fs = glob('./*/mcmc_for/03_mcmctree.out')
for infile in fs:
    gene_name = infile.split('/')[1]
    ofile = f'./pack_result/test_r2/{gene_name}.png'
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    r2 = draw_picture(infile,ofile,return_fig=False)
    collect_r2[gene_name] = r2
    
infile = '/home-user/thliao/data/nitrification_for/dating_for/rough_dating/01_mcmctree_out/out'
ofile = './pack_result/test_r2/rough_dating.png'
r2 = draw_picture(infile,ofile,return_fig=False)
collect_r2['rough dating'] = r2

infile = '/home-user/thliao/data/nitrification_for/dating_for/rough_dating/01_mcmctree_out/out'
ofile = './pack_result/test_r2/rough_dating.png'
r2 = draw_picture(infile,ofile,return_fig=False)
collect_r2['rough dating'] = r2