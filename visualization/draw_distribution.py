from os.path import *

import plotly.figure_factory as ff
from for_software.for_bayestraits.toolkit.get_result import *

df = get_df('/mnt/ivy/thliao/project/NOB/ACE/wol_phylogeny/complex_m/bst_complex.Log.txt')
xs = [list(df['q10'])]
ys = ['q10']


fig = ff.create_distplot(xs,
                            ys,
                            show_hist=True,
                            show_rug=False,
                            bin_size=0.05,
                            )

fig.layout.showlegend = False
fig.layout.height = 900
fig.layout.width = 1000
fig.layout.font.size = 25
fig.write_html(f'/mnt/ivy/thliao/project/NOB/ACE/wol_phylogeny/complex_m/q10_rate.dis.html', include_plotlyjs='cdn')

