import pandas as pd
from plotly.subplots import make_subplots

seasons = {'w':(12,1,2),
           's':(3,4,5),
           'sum':(6,7,8),
           'a': (9,10,11)}
s2color = {'w':'#2962FF',
           's':'#76FF03',
           'sum':'#D81B60',
           'a': '#FF6D00'}
#df = pd.read_csv('./大屿山marine_water_quality_north western.csv')
df = pd.read_csv('./南丫岛marine_water_quality_Southern.csv')

for rid,row in df.loc[:,['Nitrate Nitrogen (mg/L)','Nitrite Nitrogen (mg/L)','Ammonia Nitrogen (mg/L)']].iterrows():
    for col,v in row.items():
        try:
            v = float(v)
        except:
            v = float(v.replace('<',''))
            df.loc[rid,col] = v

df.loc[:,'Dates'] = pd.to_datetime(df['Dates'])
        
        
import plotly.express as px
import plotly.graph_objects as go

ndf = df.copy()
ndf = pd.DataFrame(ndf.set_index('Dates').resample('M')['DIN'].median())
ndf = ndf.loc['2009-01-31':'2020-12-31',:]
ndf.loc[:,'Dates'] = ndf.index
ndf.loc[:,'Temp'] = df.copy().set_index('Dates').resample('M')['Temperature (°C)'].median()
ndf['year'] = pd.DatetimeIndex(ndf['Dates']).year
ndf['month'] = pd.DatetimeIndex(ndf['Dates']).month

full_fig = make_subplots(2,1,shared_xaxes=True,
                        specs=[[{"secondary_y": True}], [{"secondary_y": True}]])
fig = px.scatter(ndf,x='Dates',y='DIN')
d = fig.data[0]
d['mode'] = 'markers+lines'
d['marker']['color']='#111111'
for i in [1,2]:
    full_fig.add_trace(d,i,1,secondary_y=False)

fig = px.scatter(ndf,x='Dates',y='Temp')
d = fig.data[0]
d['mode'] = 'markers+lines'
d['marker']['color']='#d23a1f'
for i in [1,2]:
    full_fig.add_trace(d,i,1,secondary_y=True)
    
####
fig = go.Figure()
max_v,min_v = 1,0
for sea,vals in seasons.items():
    ssdf = ndf.loc[ndf['month'].isin(vals),:]
    #fig.add_traces(go.Scatter(x=ssdf['Dates'],y=ssdf['DIN'],mode='markers'))
    for year,each_df in ssdf.groupby('year'):
        if sea =='w':
            half = pd.Timedelta(days=15)
            s,e = each_df.index[0],each_df.index[-2]
            s=s-half;e = e+half
            fig.add_traces( go.Scatter(x=[s,s,e,e,s], y=[min_v,max_v,max_v,min_v,min_v], fill="toself",opacity=0.7,marker=dict(color=s2color[sea]),line=dict(width=0),mode='lines',showlegend=False,hoverinfo='skip') )
            #one = pd.Timedelta(days=30)
            mid = each_df.index[-1]
            s = mid -half
            e=mid+half
            fig.add_traces( go.Scatter(x=[s,s,e,e,s], y=[min_v,max_v,max_v,min_v,min_v], fill="toself",opacity=0.7,marker=dict(color=s2color[sea]),line=dict(width=0),mode='lines',showlegend=False,hoverinfo='skip') )
        else:
            half = pd.Timedelta(days=16)
            s,e = each_df.index[0],each_df.index[-1]
            s=s-half;e = e+half
            fig.add_traces( go.Scatter(x=[s,s,e,e,s], y=[min_v,max_v,max_v,min_v,min_v], fill="toself",opacity=0.7,marker=dict(color=s2color[sea]),line=dict(width=0),mode='lines',showlegend=False,hoverinfo='skip') )
for _ in fig.data:
    full_fig.add_trace(_,1,1)


fig = go.Figure()
ssdf = ndf.loc[(ndf['month']>=3) & (ndf['month']<10),:]
#fig.add_traces(go.Scatter(x=ssdf['Dates'],y=ssdf['DIN'],mode='markers'))
for year,each_df in ssdf.groupby('year'):
    s,e = each_df.index[0],each_df.index[-1]
    fig.add_traces( go.Scatter(x=[s,s,e,e,s], y=[min_v,max_v,max_v,min_v,min_v], 
                               fill="toself",
                               opacity=0.7,
                               marker=dict(color='#7180ac'),
                               line=dict(width=0),mode='lines',showlegend=False,hoverinfo='skip') )
for _ in fig.data:
    full_fig.add_trace(_,2,1)
    
full_fig.layout.template = 'simple_white'
full_fig.layout.height = 500
full_fig.layout.width = 900
# full_fig.layout.yaxis1.title.text = 'dissolved inorgnaic <Br> nitrogen levels (mg/L)'
# full_fig.layout.yaxis2.title.text = 'dissolved inorgnaic <Br> nitrogen levels (mg/L)'
# full_fig.layout.yaxis3.title.text = 'dissolved inorgnaic <Br> nitrogen levels (mg/L)'

full_fig.update_yaxes(title_text='Dissolved inorgnaic <Br> nitrogen (DIN) (mg/L)', secondary_y=False)
full_fig.update_yaxes(title_text='Temperature (°C)', secondary_y=True)
#full_fig.layout.yaxis.title.font.color = '#7180ac'
full_fig.layout.xaxis1.range = ['2009-01-31','2020-12-31']
#full_fig.layout.xaxis.dtick='M6'
#full_fig.show()
full_fig.write_image('./南丫岛.pdf')