from collections import defaultdict

import click
import os
import plotly
import plotly.graph_objs as go
from tqdm import tqdm
from ete3 import Tree
from api_tools import tree_vis, read_tree

def get_preferred_scale(leaves1, leaves2):
    num_1 = len(leaves1)
    num_2 = len(leaves2)

    yscale = num_1 / num_2
    return 1 / yscale


def cooridinate_two_identical_trees(newick1,newick2,l2r=None):
    def sort_way(n):
        if n.is_leaf():
            return n.name
        else:
            return ';'.join(n.get_leaf_names())
    t1 = Tree(newick1,3)
    t2 = Tree(newick2,3)
    for n in t1.traverse('postorder'):
        if not n.is_leaf():
            ls = n.get_leaf_names()
            n2 = t2.get_common_ancestor([ls[0],ls[-1]])
            if set(n.get_leaf_names()) == set(n2.get_leaf_names()):
                if tuple(n.get_leaf_names()) != tuple(n2.get_leaf_names()):
                    n.children = sorted(n.children,key=lambda x:sort_way(x))
                    n2.children = sorted(n2.children,key=lambda x:sort_way(x))
    return t1,t2

def parse_color_scheme_files(file, extra_set=False, get_raw_name=False):
    """

    :param file:
    :param extra_set:
    :param get_raw_name:
    :return:
    """
    lines = open(file).read().split("\n")
    lines = [_ for _ in lines if _]
    field_labels = [_ for _ in lines if _.startswith("FIELD_LABELS")]
    field_colors = [_ for _ in lines if _.startswith("FIELD_COLORS")]
    name2color = {}
    sep_indicator = [_ for _ in lines if _.startswith("SEPARATOR")][0]
    sep_indicator = sep_indicator.split(" ")[-1]
    s2s = {"TAB": "\t", "COMMA": ",", "SPACE": " "}
    sep = s2s.get(sep_indicator, ",")

    if field_labels and field_colors:
        colors = field_colors[0].split("\t")[1:]
        anno_names = field_labels[0].split("\t")[1:]
        _name2data = [
            _ for _ in lines[lines.index("DATA") + 1 :] if _ and not _.startswith("#")
        ]
        for line in _name2data:
            line.strip("\n")
            vs = line.split("\t")
            name = vs[0]
            color = [c for c, _ in zip(colors, vs[1:]) if str(_) != "0"]
            if color:
                name2color[name] = color[0]  # would overlap.. be careful..

    else:
        lines = [
            _ for _ in lines[lines.index("DATA") + 1 :] if _ and not _.startswith("#")
        ]
        c2name = {}
        for line in lines:
            values = line.split(sep)
            name = values[0]
            if "range" in line:
                color = values[2]
                c2name[color] = values[-1]
            else:
                color = values[1]
            name2color[name] = color
        if get_raw_name:
            return c2name
    if str(extra_set) == "rename":
        new_name2color = {
            k.split("_")[-1].replace(".", "v"): v for k, v in name2color.items()
        }
        name2color = new_name2color.copy()
    return name2color


def tanglegram(
    newick1,
    newick2,
    l_color_setting=None,
    r_color_setting=None,
    l_legnth="max",
    sep="_",
    extra_set=False,
    identical=False,
    relationship_l2r=None,
):
    """

    If you want to annotate some lines descending to the particular clade
    
    For example:
    datas2 = [_ for _ in fig.data if _['xaxis'] == 'x3']  # the right tree
    for idx in right_tre.get_index('GCF_900248165.1','GCF_900696045.1'):
        datas2[idx]['line']['color']='#ff0000'  
    datas1 = [_ for _ in fig.data if _['xaxis'] == 'x']  # the left tree
    for idx in right_tre.get_index('GCF_900248165.1','GCF_900696045.1'):
        datas1[idx]['line']['color']='#ff0000'          
    :param newick1: gene tree
    :param newick2: species tree
    :param l_color_setting: itol annotating file or a dict
    :param r_color_setting: itol annotating file or a dict
    :param l_legnth:
    :param sep:
    :param extra_set:
    :param identical:
    :param relationship_l2r:
    :return:
    """
    left_tre = tree_vis(newick1, branch_length=True)
    right_tre = tree_vis(newick2, branch_length=True)

    left_leaves = [_.name for _ in left_tre.tree_obj.get_terminals()]
    right_leaves = [_.name for _ in right_tre.tree_obj.get_terminals()]
    yscale = get_preferred_scale(left_leaves, right_leaves)
    fig = plotly.tools.make_subplots(
        rows=1, cols=3, shared_yaxes=True, horizontal_spacing=0.05 / 3
    )
    datas1 = left_tre.get_plotly_data(yscale=yscale)
    datas2 = right_tre.get_plotly_data()

    # get dendrogram parts
    tqdm.write("drawing dendrograms")

    # add colors or something else into above datas
    l2color = {_[2]: "#000000" for _ in left_tre.labels_info}
    r2color = {_[2]: "#000000" for _ in right_tre.labels_info}
    if l_color_setting is not None and type(l_color_setting)==str:
        l2color.update(parse_color_scheme_files(l_color_setting, extra_set=False))
    elif l_color_setting is not None and type(l_color_setting)==dict:
        l2color.update(l_color_setting)
    if r_color_setting is not None and type(r_color_setting)==str:
        r2color.update(parse_color_scheme_files(r_color_setting, extra_set=extra_set))
    elif r_color_setting is not None and type(r_color_setting)==dict:
        l2color.update(r_color_setting)        
    # add color into generated trace
    tqdm.write("adding color")
    name2trace = {_["name"]: _ for _ in datas1 if _["name"] is not None}
    for l, color in l2color.items():
        if color != "#00000" and l in [_[2] for _ in left_tre.labels_info]:
            trace = name2trace[l]
            trace["line"]["color"] = color
    name2trace = {_["name"]: _ for _ in datas2 if _["name"] is not None}
    for r, color in r2color.items():
        if color != "#00000" and r in [_[2] for _ in right_tre.labels_info]:
            trace = name2trace[r]
            trace["line"]["color"] = color

    # add dendrogram parts into figure.
    # data1/newick1 would be the left, data2/newick2 would be the right and the leaves of it will point to left
    # for _ in datas:
    fig.add_traces(datas1, rows=[1] * len(datas1), cols=[1] * len(datas1))
    # for _ in datas2:
    fig.add_traces(datas2, rows=[1] * len(datas2), cols=[3] * len(datas2))

    # init data of middle part
    # get the y-coordinate information from below. put them into two dict.
    left_data = {
        label: y for x, y, label in left_tre.labels_info if label in left_leaves
    }
    right_data = {
        label: y for x, y, label in right_tre.labels_info if label in right_leaves
    }

    # get mapping relationship, default is from left to the right.. one to multi
    # so, the leaf names from the left tree should be the part of right, separate with `underline` or `space`
    if relationship_l2r is None:
        l2r = defaultdict(list)
        for leaf1 in left_data.keys():
            if identical:
                leaf2s = [r for r in right_data.keys() if leaf1 == r]
            else:
                leaf2s = [r for r in right_data.keys() if leaf1.split(sep)[0] == r]
            l2r[leaf1] = leaf2s
    else:
        l2r = relationship_l2r.copy()

    # init the data from above mapping dict
    c2data = defaultdict(lambda: ([], []))
    for l, color in l2color.items():
        _xs = c2data[color][0]
        _ys = c2data[color][1]
        rs = l2r.get(l, [])
        l_y = left_data.get(l, 0)
        for r in rs:
            if r not in right_data:
                continue
            r_y = right_data[r]
            _xs += [0, 1, None]
            _ys += [l_y, r_y, None]

    for color, (_xs, _ys) in c2data.items():
        trace = go.Scatter(
            x=_xs,
            y=_ys,
            mode="lines",
            line=dict(
                color=color,
            ),
            hoverinfo="none",
            showlegend=True,
        )
        fig.add_traces([trace], rows=[1], cols=[2])

    d = dict(showticklabels=False, zeroline=False)
    fig.layout.update(
        xaxis=d,
        xaxis2=d,
        xaxis3=d,
        yaxis=d,
    )
    fig.layout.xaxis3.autorange = "reversed"
    return fig,left_tre,right_tre

