import plotly.graph_objects as go
def format_fig(fig,schemes):
    _1,_2 = schemes[0],schemes[-1]
    fig.append_trace(go.Scatter(x=[_1,_1,_2,_2,_1], 
                               y=[23.2,25,25,23.2,23.2], 
                               mode='lines',
                            marker=dict(color='#bbbcbc'),
                            opacity=0.4,
                                showlegend=False,
                               fill="toself"),1,1)

    fig.update_layout(width=1300,
                      height=500,
                      yaxis1={'domain': [0, 0.75],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': True,
                                      'zeroline': True,
                              "tickfont":{"family":'Arial','size':20},
                             #         'ticks':""

                             },
                      yaxis2={'domain': [0.76, 1],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                              'showticklabels':False,
                                      'zeroline': False,
                              "range":[0,0.55],
                                      'ticks':""
                             },
                      xaxis1={'showticklabels':True,
                              'ticktext': [_.split('_')[-1] for _ in schemes] ,
                              'tickvals': schemes,
                              'tickmode':"array",
                             # "zeroline":True,
                             "tickfont":{"family":'Arial','size':20},
                             },
                     xaxis2={'showticklabels':False,'ticks':"",'showline':False},
                     template='simple_white')
    return fig