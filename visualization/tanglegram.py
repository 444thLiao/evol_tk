
from collections import defaultdict

import click
import os
import plotly
import plotly.graph_objs as go
from tqdm import tqdm

from api_tools.for_tree.format_tree import read_tree
from api_tools.for_tree.vis import main as get_plotly_data_from_newick

def get_leafs(newick):
    t = read_tree(newick, format='auto')
    return t.get_leaf_names()

def get_preferred_scale(newick1, newick2):
    t1 = read_tree(newick1, format='auto')
    t2 = read_tree(newick2, format='auto')
    num_1 = len(t1.get_leaf_names())
    num_2 = len(t2.get_leaf_names())

    yscale = num_1 / num_2
    return 1 / yscale


def parse_color_scheme_files(file, extra_set=False,get_raw_name=False):
    """

    :param file:
    :param extra_set:
    :param get_raw_name:
    :return:
    """
    lines = open(file).read().split('\n')
    lines = [_ for _ in lines if _]
    field_labels = [_ for _ in lines if _.startswith("FIELD_LABELS")]
    field_colors = [_ for _ in lines if _.startswith("FIELD_COLORS")]
    name2color = {}
    sep_indicator = [_ for _ in lines if _.startswith("SEPARATOR")][0]
    sep_indicator = sep_indicator.split(' ')[-1]
    s2s = {"TAB": '\t', "COMMA": ',', "SPACE": ' '}
    sep = s2s.get(sep_indicator, ',')

    if field_labels and field_colors:
        colors = field_colors[0].split('\t')[1:]
        anno_names = field_labels[0].split('\t')[1:]
        _name2data = [_ for _ in lines[lines.index("DATA") + 1:] if _ and not _.startswith('#')]
        for line in _name2data:
            line.strip('\n')
            vs = line.split('\t')
            name = vs[0]
            color = [c for c, _ in zip(colors, vs[1:]) if str(_) != '0']
            if color:
                name2color[name] = color[0]  # would overlap.. be careful..

    else:
        lines = [_ for _ in lines[lines.index("DATA") + 1:] if _ and not _.startswith('#')]
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
    if str(extra_set) == 'rename':
        new_name2color = {k.split('_')[-1].replace('.', 'v'): v
                          for k, v in name2color.items()}
        name2color = new_name2color.copy()
    return name2color


def main(newick1, newick2,
         color_file1=None, color_file2=None,
         l_legnth='max', sep='_', extra_set=False,
         identical=False, relationship_l2r=None):
    """

    :param newick1: gene tree
    :param newick2: species tree
    :param color_file1:
    :param color_file2:
    :param l_legnth:
    :param sep:
    :param extra_set:
    :param identical:
    :param relationship_l2r:
    :return:
    """
    left_leaves = get_leafs(newick1)
    right_leaves = get_leafs(newick2)
    yscale = get_preferred_scale(newick1, newick2)
    fig = plotly.tools.make_subplots(rows=1, cols=3, shared_yaxes=True,
                                        horizontal_spacing=0.05/3)

    # get dendrogram parts
    tqdm.write('drawing dendrograms')
    datas, labels, _, labels_draw_text, labels_x, labels_y = get_plotly_data_from_newick(newick1,
                                                                                         fixed_length=l_legnth,
                                                                                         yscale=yscale)
    datas2, labels2, _, labels_draw_text2, labels_x2, labels_y2 = get_plotly_data_from_newick(newick2,
                                                                                              fixed_length=l_legnth)
    # add colors or something else into above datas
    l2color = {_: '#000000' for _ in labels_draw_text}
    r2color = {_: '#000000' for _ in labels_draw_text2}

    if color_file1 is not None:
        l2color.update(parse_color_scheme_files(color_file1, extra_set=False))
    if color_file2 is not None:
        r2color.update(parse_color_scheme_files(color_file2, extra_set=extra_set))
    # add color into generated trace
    tqdm.write('adding color')
    name2trace = {_['name']: _ for _ in datas if _['name'] is not None}
    for l, color in l2color.items():
        if color != '#00000' and l in labels_draw_text:
            trace = name2trace[l]
            trace['line']['color'] = color
    name2trace = {_['name']: _ for _ in datas2 if _['name'] is not None}
    for r, color in r2color.items():
        if color != '#00000' and r in labels_draw_text2:
            trace = name2trace[r]
            trace['line']['color'] = color

    # add dendrogram parts into figure.
    # data1/newick1 would be the left, data2/newick2 would be the right and the leaves of it will point to left
    # for _ in datas:
    fig.add_traces(datas, rows=[1] * len(datas), cols=[1] * len(datas))
    # for _ in datas2:
    fig.add_traces(datas2, rows=[1] * len(datas2), cols=[3] * len(datas2))

    # init data of middle part
    # get the y-coordinate information from below. put them into two dict.
    left_data = {k:v
                 for k,v in dict(zip(labels_draw_text, labels_y)).items()
                 if k in left_leaves}
    right_data = {k:v
                  for k,v in dict(zip(labels_draw_text2, labels_y2)).items()
                  if k in right_leaves}

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
        trace = go.Scatter(x=_xs,
                           y=_ys,
                           mode='lines',
                           line=dict(color=color, ),
                           hoverinfo='none',
                           showlegend=True)
        fig.add_traces([trace], rows=[1], cols=[2])
    fig.layout.xaxis3.autorange = 'reversed'
    fig.layout.xaxis.showticklabels = False
    fig.layout.xaxis2.showticklabels = False
    fig.layout.xaxis3.showticklabels = False
    fig.layout.xaxis.zeroline = False
    fig.layout.xaxis2.zeroline = False
    fig.layout.xaxis3.zeroline = False
    fig.layout.yaxis.zeroline = False
    fig.layout.yaxis.showticklabels = False

    return fig