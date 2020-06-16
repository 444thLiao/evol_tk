"""
For visualization the sorted orthologous table.
raw
"""

import pandas as pd
import plotly.graph_objects as go

# raw part
infile = './sorted_OG_renamed.csv'
OG_df = pd.read_csv(infile, sep='\t', index_col=0)


def f(x):
    if not pd.isna(x):
        return 1
    else:
        return 0


val_df = OG_df.applymap(f)

fig = go.Figure(data=go.Heatmap(z=val_df.T.values,
                                colorscale=[[0, 'rgb(247, 251, 255)'],
                                            [1, 'rgb(8, 48, 107)']]))
fig.layout.width = 3000
fig.layout.height = 1000
fig.write_html('./test.html')
