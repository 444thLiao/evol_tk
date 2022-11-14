from plotly.subplots import make_subplots
import plotly.graph_objects as go
from dating_workflow.debug_for.draw_infinite_plot import fit_line,read_multi_mcmc


## need to be adjusted
rename_f = lambda x:x
mcmc_list = ['/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/dating/plastids_test/mcmc_o/B1E6.anammox_euk_IR_prot/mcmc_for/mcmc.txt']
s2color = {"cog25": "#004D40", "mito": "#FBC02D", "MERGE": "#B71C1C"}
post_dfdict, post_CIs_dict = read_multi_mcmc(mcmc_list, rename_f=rename_f)
schemes = sorted(schemes,key=lambda x: (int(x.split(';')[1].strip('C')), s.index(x.split(';')[2])))
target_nodes = []
###### 
fig = make_subplots(rows=2, cols=1, shared_xaxes=True)
for key in schemes :
    df = post_dfdict[key]
    topo, cal, gene = key.split(";")
    for n in target_nodes:
        sdf = df[n]
        fig.add_trace(go.Box(
            y=sdf,
            name=key,
            marker={"color": s2color[gene]},
            boxpoints=False,
            showlegend=False,),1,1
        )
xs = schemes
ys = []
for key in schemes:
    coef, r2 = fit_line(
                x=post_CIs_dict[key]["Posterior mean time (100 Ma)"].values[:-1],   # remove lnL
                y=post_CIs_dict[key]["CI_width"].values[:-1]  # remove lnL
            )
    ys.append(coef)
fig.append_trace(go.Bar(x=xs,y=ys,
                        marker_color='#a87900',
                        text = ["{:.3f}".format(_) for _ in ys],
                        showlegend=False,textposition='outside'),2,1)
fig.show()