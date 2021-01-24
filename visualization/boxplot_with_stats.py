"""
Provide framework for visualizing boxplot with stats comparison

"""
import pandas as pd 
from plotly.subplots import make_subplots
from scipy.stats import mannwhitneyu
import plotly.graph_objects as go
import itertools 

def get_stars(p):
    if p <= 0.05 and p > 0.01:
        return '*'
    elif p <= 0.01 and p > 0.001:
        return '**'
    elif p <= 0.001:
        return '***'
    return 'ns'


draw_df = pd.DataFrame()

x_used_col = 'env type'
x_iter = draw_df[x_used_col].unique()

compared_col = 'nifH type'
f1_f2_sets = list(itertools.combinations(draw_df[compared_col].unique(),2))

compared_v_col = 'ratio'


fig = make_subplots(rows=1, cols=4, shared_yaxes=True, horizontal_spacing=0.01,
                    subplot_titles=['marine (n)', 'soil', 'others', 'plant']
                    )
for i, env in enumerate(['marine', 'soil', 'others', 'plant']):
    traces = []
    sub_df = draw_df.loc[draw_df['env type'] == env, :]
    for nifH_t, idxs in sub_df.groupby('nifH type').groups.items():
        _sub_df = sub_df.loc[sub_df['nifH type'] == nifH_t, :]
        traces.append(go.Violin(x=[nifH_t]*len(idxs),
                                y=_sub_df.loc[idxs, 'ratio'],
                                legendgroup=nifH_t,
                                name=nifH_t,
                                line_color=t2color[nifH_t],
                                spanmode='hard',
                                width=0.7,
                                points=False,
                                meanline=dict(visible=True),
                                showlegend=False
                                ))
    fig.add_traces(traces, 1, i+1)
    
for idx, _xi in enumerate(x_iter):
    pos = [1.5, 0.5, 1, 2, 3, 4, 5]
    # need to manual adjust
    base_height = 110
    # dependent on your data
    per_height = 5 # used to shift along with y-axis
    _count = 0  # used to shift along with x-axis
    for f1, f2 in f1_f2_sets:
        r1 = draw_df.loc[(draw_df[x_used_col] == _xi) &
                         (draw_df[compared_col] == f1), compared_v_col]
        r2 = draw_df.loc[(draw_df[x_used_col] == _xi) &
                         (draw_df[compared_col] == f2), compared_v_col]
        t = mannwhitneyu(r1, r2,
                         )
        p = t.pvalue

        fig.add_traces(go.Scatter(x=[f1, f1, f2, f2],
                                  y=[base_height,
                                     base_height+per_height,
                                     base_height+per_height,
                                     base_height],
                                  showlegend=False,
                                  mode='lines', line=dict(color='#000000')), 1, idx+1)
    #     fig.add_annotation(
    #         x=pos[_count],
    #         y=base_height+12,
    #         xanchor='center',
    #         xref=f'x{idx+1}',
    #         font=dict(size=15),
    #         showarrow=False,

    #         text="<Br> {:.2e} ".format(p),
    #     )
        fig.add_annotation(
            x=pos[_count],
            y=base_height+10,
            xanchor='center',
            xref=f'x{idx+1}',
            font=dict(size=15),
            showarrow=False,
            text=get_stars(p),
        )
        #print(get_stars(p),p,env,f1,f2)
        base_height += 10
        _count += 1
        

fig.layout.template = "simple_white"
fig.layout.width = 1400
fig.layout.height = 500
fig.layout.font.size = 20
#fig.layout.yaxis.title = "Relative abundance of <Br> different types of nif clusters"
for _ in fig.layout.annotations[:4]:
    _['font']['size'] = 23
fig.layout.yaxis.tickvals = [0, 20, 40, 60, 80, 100]
fig.layout.yaxis.title = "Relative abundance (%)"
# fig.show()