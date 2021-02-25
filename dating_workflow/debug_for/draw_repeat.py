import pandas as pd
from glob import glob
import plotly.express as px
import plotly.graph_objects as go
from os.path import *
import os
from dating_workflow.debug_for.draw_infinite_plot import get_CI

def get_r2(text):
    r_squre_text = text['hovertemplate'].split('<br>')[2]
    coef = text['hovertemplate'].split('<br>')[1]
    coef = coef.replace('time2','y')
    coef = coef.replace('time1', 'x').split('+')[0].strip()
    prefix, v = r_squre_text.split('=')

    return f"{coef}, {prefix} = {v}"


def draw_r(time1,time2):
    df = pd.DataFrame()
    df.loc[:,'time1'] = time1/10
    df.loc[:, 'time2'] = time2/10
    fig = px.scatter(df,
                     x='time1',
                     y='time2',
                     trendline="ols")

    r_squre_text = [get_r2(_)
                        for _ in [fig.data[1]]]

    fig.update_layout(
        showlegend=False,
        annotations=[
            go.layout.Annotation(
                x=0.5,
                y=max(list(df.loc[:, 'time1'])),
                text='<br>'.join(r_squre_text),
                showarrow=False,
                font=dict(size=15)
            )
        ]
    )
    fig.layout.xaxis.title.text = 'Posterior mean time (Gya) − Run1'
    fig.layout.yaxis.title.text = 'Posterior mean time (Gya) − Run2'
    return fig


def main(pattern,odir):
    if not exists(odir):
        os.makedirs(odir)
    a = glob(pattern)
    for f1 in a:
        if '187' not in f1:
            name = f1.split('_')[2]
            set_name = f1.split('_')[3]
            f2 = f1.replace('_run1','_run2')
            df1,df2 = get_CI(f1),get_CI(f2)
            time1,time2 = df1.iloc[:,-1],df2.iloc[:,-1]
            fig = draw_r(time1,time2)
            fig.write_image(join(odir,f'repeat_{name}_{set_name}.png'))

if __name__ == '__main__':
    odir = './dating_for/result_draw'
    pattern = "./dating_for/clock3/*_run1/run.log"
    main(pattern,odir)