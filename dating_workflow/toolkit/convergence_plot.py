import plotly.express as px
import pandas as pd
import plotly.graph_objs as go

def parse_mcmc(infile):
    t = pd.read_csv(infile,sep='\t',index_col=0)
    return t.mean()

def get_r2(text):
    r_squre_text = text['hovertemplate'].split('<br>')[2]
    prefix, v = r_squre_text.split('=')
    # mode_text = text['hovertemplate'].split('<br>')[4]
    # mode = mode_text.split('=')[-1]

    return f"{prefix} ={v}"


mc1 = './dating_for/198g_3cal_set3_rgene_1_10/mcmc_for/mcmc.txt'
mc2 = './dating_for/198g_3cal_set3_rgene_1_10/mcmc_for/repeat04/mcmc.txt'
def draw_conv_plot(mc1,mc2):
    t1 = parse_mcmc(mc1)
    t2 = parse_mcmc(mc2)
    nodes = [_ for _ in t1.index if _.startswith('t_')]
    _df = pd.DataFrame({'run 01':t1,'run 02':t2})
    _df = _df.loc[nodes,:]
    fig = px.scatter(_df,x='run 01',y='run 02',trendline="ols")
    r_squre_text = [get_r2(fig.data[1])]
    fig.update_layout(
        showlegend=False,
        annotations=[
            go.layout.Annotation(
                x=0,
                y=max(list(_df.loc[:, 'run 01'])),
                text='<br>'.join(r_squre_text),
                showarrow=False,
                font=dict(size=20)
            )
        ]
    )
    fig.layout.width = 1000
    fig.layout.height = 1000
    fig.write_html('./test.html')