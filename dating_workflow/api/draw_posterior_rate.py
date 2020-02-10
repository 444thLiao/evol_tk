import plotly.express as px

# from .IO_mcmctree import read_outfile
f = "./dating_for/83g/83g_set1/mcmc_for/03_mcmctree.out"
df = read_outfile(f)
df = df.reindex([_ for _ in df.index if _.startswith('r_g')])
df.columns = ['CI_width', 'CIs', 'Posterior mean rate (per AA per 100 Mya)']
df.loc[:,'gene'] = [_.split("_")[1].strip('g')
                    for _ in df.index]
df.loc[:,'num_branch'] = [_.split('_')[2]
                          for _ in df.index]
fig = px.scatter(df,x='num_branch',
                 y='Posterior mean rate (per AA per 100 Mya)',
                 color="gene",
                 marginal_y="violin")
fig.write_html('./test.html')

fig.write_image('./test.png')