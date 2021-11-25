import plotly
import plotly.graph_objects as go
from tanglegram import tree_vis

gid2l = {
    row.split(",")[0]: row.split(",")[2]
    for row in open("/home-user/mxie/wgs_test/16S/location.txt").read().split("\n")
    if len(row.split(",")) == 3
}
gid2c = {
    row.split(",")[0]: row.split(",")[2]
    for row in open("/home-user/mxie/wgs_test/16S/compartment.txt").read().split("\n")
    if len(row.split(",")) == 3
}
gid2g = {
    row.split(",")[0].split("_")[0]: row.split(",")[2]
    for row in open("/home-user/mxie/wgs_test/16S/group-16S.txt").read().split("\n")
    if len(row.split(",")) == 3
}

cmap_g = {"sister group": "#689F38", "target": "#0277BD", "outgroup": "#E65100"}
cmap_c = {"mucus": "#D3D3D3", "skeleton": "#A9A9A9"}
cmap_l = {
    "Bluff Island": "#E65100",
    "Lo Chau": "#689F38",
    "Sham Wan": "#0277BD",
    "Yam Tsai Wan": "#000000",
}


def get_tre(infile):
    tre = Tree(infile)
    gids = tre.get_leaf_names()
    dropped_leaves = []
    for n in "L3,P4,O8,I2,C3".split(","):
        if n in gids:
            dropped_leaves.append(n)
    final_gids = set(gids).difference(set(dropped_leaves))
    tre.prune(final_gids)

    rooted_lca = "GNM003217735|GNM007923355".split("|")
    tre.set_outgroup(tre.get_common_ancestor(rooted_lca))
    tre = sort_tree(tre)
    ordered_gids = tre.get_leaf_names()
    return tre


from collections import defaultdict
from tqdm import tqdm
import pandas as pd

# g2g_mash_v = defaultdict(dict)
# for row in tqdm(open("/home-user/mxie/wgs_test/data/pairwise_merged_db.dist")):
#     s1, s2, mash_v, p, hashed = row.split('\t')
#     g1 = s1.split('/')[-1].replace('.fna', '')
#     g2 = s2.split('/')[-1].replace('.fna', '')
#     g2g_mash_v[g1][g2] = float(mash_v)
# df = pd.DataFrame.from_dict(g2g_mash_v)
# df.to_csv(f"/mnt/ivy/thliao/project/coral_ruegeria/data_processing/pairwise_merged_db.tab", sep='\t', index=1)

tre = get_tre(
    "/mnt/ivy/thliao/project/coral_ruegeria/data_processing/phylogeny/iqtree/bac120_over20p_trimal.contree"
)
df = pd.read_csv(
    "/mnt/ivy/thliao/project/coral_ruegeria/data_processing/pairwise_merged_db.tab",
    sep="\t",
    index_col=0,
)
sub_gids = [_ for _ in ordered_gids if _ in df.index]
tre.prune(sub_gids)
tre = sort_tree(tre)
ordered_gids = tre.get_leaf_names()
df = df.reindex(ordered_gids, columns=ordered_gids)

vt = tree_vis(tre, leaves2top=False)
fig = plotly.tools.make_subplots(
    rows=1, cols=3, shared_yaxes=True, horizontal_spacing=0.05 / 3
)

for gid in df.index:
    fig.add_trace(
        go.Bar(
            x=[1],
            y=[1],
            marker={
                "color": cmap_c.get(gid2c.get(gid, ""), "#ffffff"),
                "line": {"width": 0},
            },
            name=gid,
            showlegend=False,
        ),
        1,
        2,
    )

y = [vt.clade_y_pos[vt.get_clade(c)] for c in df.columns]
heatmap = go.Heatmap(
    x=df.index,
    y=y,
    z=1 - df.values,
    colorscale="YlGnBu",
    reversescale=True,
    showscale=False,
)
fig.append_trace(heatmap, 1, 3)
datas1 = vt.get_plotly_data(yscale=1, fix_length=sorted(vt.clade_x_pos.values())[-1])
fig.add_traces(datas1, rows=[1] * len(datas1), cols=[1] * len(datas1))

fig.update_layout(
    barmode="stack",
    width=1500,
    height=1500,
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    xaxis={
        "domain": [0, 0.1],
        "showgrid": False,
        "zeroline": False,
        "showticklabels": False,
        "showspikes": False,
    },
    xaxis2={
        "domain": [0.1, 0.12],
        "zeroline": False,
        "showticklabels": False,
        "showspikes": False,
    },
    xaxis3={
        "domain": [0.12, 1],
        "zeroline": False,
    },
    yaxis={
        "zeroline": False,
        "showgrid": False,
        "showticklabels": False,
        "showspikes": False,
    },
    yaxis2={
        "zeroline": False,
        "showgrid": False,
        "showticklabels": False,
        "showspikes": False,
    },
    yaxis3={
        "zeroline": False,
        "showgrid": False,
        "side": "right",
        "ticktext": df.columns,
        "tickvals": y,
        "showticklabels": True,
        "showspikes": True,
        "range": [0, sorted(vt.clade_y_pos.values())[-1]],
    },
)
fig.show()