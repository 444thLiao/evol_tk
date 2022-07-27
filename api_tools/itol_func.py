import os
from ete3 import Tree
import sys
from os.path import join, exists, dirname

import matplotlib as mpl
import numpy as np
import plotly.express as px

sys.path.insert(0, dirname(__file__))

home = os.getenv("HOME")
indir = f"{home}/template_txt/"
if not exists(indir):
    indir = join(dirname(__file__), "itol_template")

color_style_template = join(indir, "colors_styles_template.txt")
color_strip_template = join(indir, "dataset_color_strip_template.txt")
dataset_styles_template = join(indir, "dataset_styles_template.txt")
dataset_binary_template = join(indir, "dataset_binary_template.txt")
label_template = join(indir, "labels_template.txt")
dataset_symbol_template = join(indir, "dataset_symbols_template.txt")
matrix_like_template = join(indir, "dataset_external_shapes_template.txt")
dataset_text_template = join(indir, "dataset_text_template.txt")
dataset_gradient_template = join(indir, "dataset_gradient_template.txt")
dataset_piechart_template = join(indir, "dataset_piechart_template.txt")
dataset_simplebar_template = join(indir, "dataset_simplebar_template.txt")
collapse_template = join(indir, "collapse.txt")
dataset_connections_template = join(indir, "dataset_connections_template.txt")



def get_used_sep(text):
    separator = [_ for _ in text.split("\n") if _.startswith("SEPARATOR")]
    assert len(separator) == 1
    separator = separator[0].strip()
    sep = separator.split(" ")[1]
    if sep == "TAB":
        return "\t"
    elif sep == "SPACE":
        return " "
    elif sep == "COMMA":
        return ","
    
def replacing_params(text, kwarg={}):
    sep= get_used_sep(text)
    for k, v in kwarg.items():
        if k.upper() in ['MARGIN',"STRIP_WIDTH"]:
            row = [_ for _ in text.split('\n') if k.upper() in _]
            assert len(row) == 1
            if row:
                text = text.replace(row[0],sep.join([k.upper(),str(v)])+'\n' )
        else:
            text = text.replace(k, v)
    return text

def deduced_legend(info2color, info_name="dataset", sep=","):
    # for implemented a legend with dictinonary named info2color.
    legend_title = info_name
    legend_shape = sep.join(["1"] * len(info2color))
    legend_colors = sep.join([_ for _ in info2color.values()])
    legend_labels = sep.join(list(info2color.keys()))
    legend_text = f"""
LEGEND_TITLE{sep}{legend_title}
LEGEND_SHAPES{sep}{legend_shape}
LEGEND_COLORS{sep}{legend_colors}
LEGEND_LABELS{sep}{legend_labels}"""
    return legend_text


def deduced_legend2(info2style, infos, same_colors=False, sep="\t"):
    # for info2style instead of info2color
    
    colors_theme = px.colors.qualitative.Dark24 + px.colors.qualitative.Light24
    shapes = []
    labels = []
    colors = []
    for idx, info in enumerate(infos):
        shapes.append(info2style[info].get("shape", "1"))
        labels.append(info2style[info].get("info", info))
        if not same_colors:
            colors.append(info2style[info].get("color", colors_theme[idx]))
        else:
            colors.append(info2style[info].get("color", same_colors))
    legend_text = [
        "FIELD_SHAPES" + sep + sep.join(shapes),
        "FIELD_LABELS" + sep + sep.join(labels),
        "FIELD_COLORS" + sep + sep.join(colors),
    ]
    legend_text = "\n".join(legend_text)
    return legend_text


def annotate_outgroup(
    ID2infos,
    info2style,
):
    for ID, infos in ID2infos.items():
        # arbitary to choose a column of infos as template for deduced legend
        break
    template_text = deduced_legend2(info2style, infos)
    annotate_text = ""
    for ID, infos in ID2infos.items():
        row = [ID]
        for info in infos:
            if info not in info2style:
                row.append("-1")
            elif info in info2style:
                row.append(info2style[info]["status"])
        annotate_text += "\t".join(row) + "\n"
    return template_text + annotate_text


def to_simple_bar(id2val, dataset_name="Value", other_params={}):
    # id2value, could be {ID: number /string like number}
    used_template = dataset_simplebar_template
    template_text = open(used_template).read().format(
        dataset_name=dataset_name) + "\n"
    sep = get_used_sep(template_text)
    c = []
    for gid, val in id2val.items():
        val = str(val)
        c.append(f"{gid}{sep}{val}")
    template_text = replacing_params(template_text, other_params)
    return template_text + "\n".join(c)

def to_connections(a2b_tuple_list,pair2style={},
                   default_color="#007acc",
                   default_width="1",
                   default_style="normal",
                   default_label="",
                   dataset_name="HGT transfers",
                   other_params={}):
    used_template = dataset_connections_template
    template_text = open(used_template).read() + "\n"
    sep = get_used_sep(template_text)
    annotate_texts = []
    for (a,b) in a2b_tuple_list:
        pair = (a,b)
        width = pair2style.get(pair,{}).get('width',default_width)
        color = pair2style.get(pair,{}).get('color',default_color)
        style = pair2style.get(pair,{}).get('style',default_style)
        label = pair2style.get(pair,{}).get('label',default_label)
        annotate_texts.append(sep.join([a,b,width,color,style,label]))
    annotate_text = '\n'.join(annotate_texts)
    
    template_text = template_text.format(
        dataset_label=dataset_name, legend_text="", 
    )
    template_text = replacing_params(template_text, other_params)
    return template_text + annotate_text

def to_binary_shape(
    ID2info,
    info2style=None,
    same_color=False,
    dataset_name="Presence/Absence matrix",
    manual_v=[],
    unfilled_other=False,
    other_params={},
    no_legend=False,
    one_row=False,
):
    # id2info, could be {ID:list/set}
    # info2color: could be {gene1: {shape:square,color:blabla},}
    # None will use default.
    # if turn unfilled_other on, it will not draw the unfilled markers
    #
    template_text = open(dataset_binary_template).read() + "\n"
    sep = get_used_sep(template_text)
    if not manual_v:
        all_v = list(sorted(set([_ for v in ID2info.values() for _ in v if _])))
    else:
        all_v = manual_v
    if info2style is None:
        info2style = {k: {} for k in all_v}
    if one_row and len(set([_ for v in ID2info.values() for _ in v])) ==2: # give one_row params and only one value
        unfilled_label = '0'
    else:
        unfilled_label = "-1" if unfilled_other else "0"

    annotate_text = []
    for ID, vset in ID2info.items():
        row = sep.join(
            [ID] + ["1" if _ in vset else unfilled_label for _ in all_v])
        annotate_text.append(row)
    annotate_text = "\n".join(annotate_text)

    if no_legend:
        real_legend_text = ""
    else:
        legend_text = deduced_legend2(
            info2style, all_v, sep=sep, same_colors=same_color)
        real_legend_text = f"LEGEND_TITLE\t{dataset_name}\n" + legend_text.replace("FIELD", "LEGEND")
    template_text = replacing_params(template_text, other_params)
    template_text = template_text.format(
        legend_text=legend_text + "\n" + real_legend_text, dataset_label=dataset_name
    )
    return template_text + "\n" + annotate_text


def to_color_strip(ID2info, info2color, dataset_name="strip", other_params={}):
    # id2info, could be {ID: name}
    # info2color: could be {name: color,}
    ID2info = {k: str(v) for k, v in ID2info.items()}

    template_text = open(color_strip_template).read()
    sep = get_used_sep(template_text)

    id2color = {id: info2color[info] for id, info in ID2info.items()}
    annotate_text = "\n".join(
        [sep.join([id, color, ID2info[id]]) for id, color in id2color.items()]
    )
    legend_text = deduced_legend(info2color, dataset_name)
    dataset_name = dataset_name.replace("/", "_")

    template_text = replacing_params(template_text, other_params)
    template_text = template_text.format(
        legend_text=legend_text, dataset_label=dataset_name
    )
    return template_text + "\n" + annotate_text


def to_color_labels_bg(ID2info, info2color, dataset_name="labels bg", other_params={}):
    # clade for
    # id2info, could be {ID: name}
    # info2color: could be {name: color,}
    template_text = open(dataset_styles_template).read()
    sep = get_used_sep(template_text)

    id2color = {ID: info2color[info] for ID, info in ID2info.items()}

    # row_remplate = sep.join(["{ID}\t{TYPE}\t{WHAT}\t{COLOR}\t{WIDTH_OR_SIZE_FACTOR}\t{STYLE}\t{BACKGROUND_COLOR}"])
    legend_text = deduced_legend(info2color, dataset_name, sep="\t")

    template_text = template_text.format(
        dataset_label=dataset_name, legend_text="")

    rows = [
        sep.join([ID, "label", "node", "#000000", "1", "bold", color])
        for ID, color in id2color.items()
    ]
    template_text = replacing_params(template_text, other_params)
    return template_text + "\n".join(rows)


def to_color_branch(ID2info, info2color, dataset_name="color branch", no_legend=False):
    # clade for
    # id2info, could be {ID: name}
    # info2color: could be {name: color,}
    template_text = open(dataset_styles_template).read()
    id2col = {ID: info2color[info] for ID, info in ID2info.items()}
    each_template = "{ID}\t{TYPE}\t{WHAT}\t{COLOR}\t{WIDTH_OR_SIZE_FACTOR}\t{STYLE}\t{BACKGROUND_COLOR}"
    if no_legend:
        legend_text = ""
    else:
        legend_text = deduced_legend(info2color, dataset_name, sep="\t")

    template_text = template_text.format(
        dataset_label=dataset_name, legend_text=legend_text
    )
    rows = [
        each_template.format(
            ID=ID,
            TYPE="branch",
            WHAT="node",
            COLOR=color,
            WIDTH_OR_SIZE_FACTOR=3,
            STYLE="normal",
            BACKGROUND_COLOR="",
        )
        for ID, color in id2col.items()
    ]
    rows += [
        each_template.format(
            ID=ID,
            TYPE="label",
            WHAT="node",
            COLOR=color,
            WIDTH_OR_SIZE_FACTOR=1,
            STYLE="bold",
            BACKGROUND_COLOR="",
        )
        for ID, color in id2col.items()
    ]
    return template_text + "\n".join(rows)


def to_color_range(
    ID2info, info2color, dataset_name="color branch", no_legend=True, other_params={}
):
    # add color range
    # id2info, could be {ID: name}
    # info2color: could be {name: color,}
    template_text = open(color_style_template).read()
    sep = get_used_sep(template_text)
    id2color = {ID: info2color[info] for ID, info in ID2info.items()}
    # each_template = "{ID}\trange\t{COLOR}\t{name}"
    if no_legend:
        legend_text = ""
    else:
        legend_text = deduced_legend(info2color, dataset_name, sep=sep) + "\n"

    template_text = template_text.format(
        dataset_label=dataset_name, legend_text=legend_text
    )
    rows = [sep.join([ID, 'range',color, str(ID2info[ID])])
            for ID, color in id2color.items()]
    template_text = replacing_params(template_text, other_params)
    return template_text + "\n".join(rows)


def to_color_Clade(
    ID2info,
    info2color,
    tree,
    dataset_name="color branch clade",
    no_legend=False,
    bgcolor={},
):
    # colorize branch within whole clade when all children under this clade share same info.
    template_text = open(dataset_styles_template).read()
    id2col = {ID: info2color[info] for ID, info in ID2info.items()}
    each_template = "{ID}\t{TYPE}\t{WHAT}\t{COLOR}\t{WIDTH_OR_SIZE_FACTOR}\t{STYLE}\t{BACKGROUND_COLOR}"
    if no_legend:
        legend_text = ""
    else:
        legend_text = deduced_legend(info2color, dataset_name, sep="\t")

    tree_obj = Tree(tree, format=3)

    def collapsed_leaf(node):
        if node.is_leaf():
            return True
        else:
            leafs_all = node.get_leaves()
            children_names = [ID2info.get(_.name, "") for _ in leafs_all]
            if "" in children_names:
                return False
            if len(set(children_names)) == 1:
                return True
            else:
                return False

    new_tree_obj = Tree(tree_obj.write(is_leaf_fn=collapsed_leaf))
    new_leaves_names = [_.name for _ in new_tree_obj.get_leaves()]
    internal_nodes = [tree_obj.search_nodes(
        name=_)[0] for _ in new_leaves_names]
    internal_nodes = [_ for _ in internal_nodes if not _.is_leaf()]
    internal_node2info = {
        n.name: info2color.get(ID2info.get(n.get_leaves()[0].name, ""), "")
        for n in internal_nodes
    }
    internal_node2info = {
        n.name: info2color[ID2info[n.get_leaves()[0].name]] for n in internal_nodes
    }
    # dropped_IDs = [_.name for node in internal_nodes for _ in node.get_leaves()]
    # for ID in list(id2col.keys()):
    #     id2col.pop(ID)

    template_text = template_text.format(
        dataset_label=dataset_name, legend_text=legend_text
    )
    rows = [
        each_template.format(
            ID=ID,
            TYPE="branch",
            WHAT="node",
            COLOR=color,
            WIDTH_OR_SIZE_FACTOR=7,
            STYLE="normal",
            BACKGROUND_COLOR="",
        )
        for ID, color in id2col.items()
    ]

    rows += [
        each_template.format(
            ID=ID,
            TYPE="label",
            WHAT="node",
            COLOR=color,
            WIDTH_OR_SIZE_FACTOR=1,
            STYLE="bold",
            BACKGROUND_COLOR=bgcolor.get(ID, ""),
        )
        for ID, color in id2col.items()
    ]
    rows += [
        each_template.format(
            ID=ID,
            TYPE="branch",
            WHAT="clade",
            COLOR=color,
            WIDTH_OR_SIZE_FACTOR=7,
            STYLE="normal",
            BACKGROUND_COLOR="",
        )
        for ID, color in internal_node2info.items()
    ]
    return template_text + "\n".join(rows)


def to_node_symbol(in_tree, dataset_name="bootstrap", maxsize="15", other_params={}):
    # normally for bootstrap
    # give it a tree is enough, must have internal node name.
    template_text = open(dataset_symbol_template).read()
    sep = get_used_sep(template_text)

    from api_tools.for_tree.format_tree import read_tree
    tree = read_tree(in_tree, format=3)
    id2node = {n.name: n for n in tree.traverse()}
    id2support = {}
    for n in tree.traverse():
        if (not n.is_leaf()) and n.name:
            try:
                support_v = int(n.name.split("_S")[-1])
            except:
                support_v = 100
            id2support[n.name] = support_v
    # ID,symbol,size,color,fill,position,label

    rows = []
    for id, s_v in id2support.items():
        size = "7"
        shape = "2"
        filled = "1"
        if int(s_v) >= 95:
            color = "#000000"
        elif int(s_v) >= 85 and int(s_v) < 95:
            color = "#777777"
        elif int(s_v) >= 65 and int(s_v) < 85:
            color = "#eeeeee"
        else:
            color = ""
        if color:
            node = id2node[id]
            new_id = (
                node.children[0].get_leaf_names()[0]
                + "|"
                + node.children[1].get_leaf_names()[0]
            )
            row = sep.join([new_id, shape, size, color, filled, "1", ""])
            rows.append(row)

    annotate_text = "\n".join(rows)
    template_text = template_text.format(
        dataset_label=dataset_name, legend_text="", maximum_size=maxsize
    )
    template_text = replacing_params(template_text, other_params)
    return template_text + annotate_text


def get_text_anno(id2val, other_params={}, dataset_name="numerical text"):
    """
    It could generate text visualization for color gradient
    Args:
        id2val ([dict]): [description]
        other_params ([dict]): [description]

    Returns:
        [type]: [description]
    """

    template_text = open(matrix_like_template).read()
    sep = get_used_sep(template_text)

    template_text = template_text.format(
        dataset_label=dataset_name,
        # must use the color with alpha/transparency
        field_color="rgba(0,255,0,0)",
        field_labels=sep.join(["text"]),
    )
    annotate_text = []
    for id, val in id2val.items():
        annotate_text.append(f"{id}{sep}{val:.2f}")

    template_text = replacing_params(template_text, other_params)
    # template_text = template_text.replace("#HEIGHT_FACTOR,1","HEIGHT_FACTOR\t1.5")
    return template_text + "\n".join(annotate_text)


def to_matrix_shape(ID2categorical_v, dataset_name, color="#000000", other_params={}):
    # id2info, could be {ID: name}

    template_text = open(matrix_like_template).read()
    sep = get_used_sep(template_text)

    all_v = set(map(str, ID2categorical_v.values()))
    all_v = list(sorted(all_v))
    if type(color) == str:
        color_str = sep.join([color] * len(all_v))
    elif type(color) == dict:
        color_str = sep.join([color.get(_, "#000000") for _ in all_v])
    else:
        raise IOError
    template_text = template_text.format(
        dataset_label=dataset_name,
        field_color=color_str,
        field_labels=sep.join(all_v),
    )
    annotate_text = ""
    for ID, v in ID2categorical_v.items():
        vals = [ID]
        for _ in range(len(all_v)):
            if v == str(all_v[_]):
                vals.append("10")
            else:
                vals.append("0")
        annotate_text += sep.join(vals) + "\n"
    template_text = replacing_params(template_text, other_params)
    return template_text + annotate_text


def to_label(id2new_id, other_params={}):
    template_text = open(label_template).read()
    sep = get_used_sep(template_text)
    full_text = replacing_params(template_text, other_params)

    full_text = template_text[::]
    for id, new_id in id2new_id.items():
        id = str(id)
        new_id = str(new_id)
        full_text += sep.join([id, new_id])+'\n'  # "%s,%s\n" % (id, new_id)
    return full_text


def colorFader(c1, c2, mix=0
               ):  # (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1 = np.array(mpl.colors.to_rgb(c1))
    c2 = np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1 - mix) * c1 + mix * c2)


# generate_gradient_legend(100,50,0,'#ff0000','#FFFFFF','#0000ff')
def generate_gradient_legend(
    max_val, mid_val, min_val, max_c, mid_c, min_c, num_interval=7
):
    legened_v2color = {}
    if num_interval % 2 == 1:
        remained_i = (num_interval - 1) // 2
        legened_v2color[round(mid_val, 2)] = mid_c
    else:
        remained_i = num_interval // 2

    total_v = max_val - mid_val
    inter_v = int(total_v / remained_i)
    for _p in range(1, remained_i):
        des_v = _p * inter_v + mid_val
        per = _p * inter_v / total_v
        new_color = colorFader(mid_c, max_c, per)
        legened_v2color[round(des_v, 2)] = new_color
    legened_v2color[round(max_val, 2)] = max_c

    total_v = mid_val - min_val
    inter_v = int(total_v / remained_i)
    for _p in range(1, remained_i):
        des_v = _p * inter_v + min_val
        per = _p * inter_v / total_v
        new_color = colorFader(min_c, mid_c, per)
        legened_v2color[round(des_v, 2)] = new_color
    legened_v2color[round(min_val, 2)] = min_c
    return legened_v2color


def color_gradient(
    id2val, dataset_name="Completness", max_val=None, min_val=None, mid_val=50, other_params={}
):
    default_max = "#ff0000"
    default_min = "#0000ff"
    default_mid = "#FFFFFF"

    template_text = open(dataset_gradient_template).read()
    sep = get_used_sep(template_text)

    all_vals = list(set([v for k, v in id2val.items()]))

    mid_val = np.mean(all_vals) if mid_val is None else mid_val
    max_val = max(all_vals) if max_val is None else max_val
    min_val = min(all_vals) if min_val is None else min_val

    l2colors = generate_gradient_legend(
        max_val, mid_val, min_val, default_max, default_mid, default_min, num_interval=7
    )

    legend_text = f"""
LEGEND_TITLE{sep}{dataset_name}
LEGEND_SHAPES{sep}{sep.join(['1'] * 7)}
LEGEND_COLORS{sep}{sep.join([_[1] for _ in list(sorted(l2colors.items()))])}
LEGEND_LABELS{sep}{sep.join(map(str, [_[0] for _ in list(sorted(l2colors.items()))]))}"""

    annotate_text = "\n".join(
        [f"{label}{sep}{val}" for label, val in id2val.items()])

    text = template_text.format(
        dataset_label=dataset_name,
        legend_text=legend_text,
        color_min=default_min,
        color_max=default_max,
        color_mid=default_mid)
    text = replacing_params(text, other_params)

    return text + "\n" + annotate_text

def pie_size_chart(id2val, color='#ff0000', 
                   dataset_name="habitat prob", pos=1, factor=1):
    template_text = open(dataset_piechart_template).read()
    sep = get_used_sep(template_text)
    annotate_text = []
    
    for gid in id2val:
        cat_vals = [gid, pos, str(int(id2val.get(gid,0)*factor))]
        cat_vals.append('1')
        cat_vals = sep.join([str(_) for _ in cat_vals])
        # gid,pos = str(gid),str(pos)
        annotate_text.append(cat_vals)

    field_labels = sep.join(['Size'])
    field_colors = sep.join([color])
    template_text = template_text.format(
        dataset_label=dataset_name,
        field_colors=field_colors,
        field_labels=field_labels,
        legend_text=''
    )
    final_text = template_text + "\n" + "\n".join(annotate_text)
    return final_text


def pie_chart(id2cat2val, cat2style, dataset_name="habitat prob",scaled_radius=False, pos=1,legend=False):
    """
    :param id2cat2val:
    :param cat2style:
    :param dataset_label:
    :return:
    """
    template_text = open(dataset_piechart_template).read()
    sep = get_used_sep(template_text)
    annotate_text = []
    all_cat = [_k for _k, v in id2cat2val.items() for _k in v]
    all_cat = list(set(all_cat))
    sorted_cat = sorted(all_cat)
    
        
    for gid in id2cat2val:
        if not scaled_radius:
            d = 10
        else:
            d = sum(id2cat2val[gid].values())
        cat_vals = [gid, pos, str(d)]
        for cat in sorted_cat:
            cat_vals.append(str(id2cat2val[gid].get(cat, "0")))
        cat_vals = sep.join([str(_) for _ in cat_vals])
        # gid,pos = str(gid),str(pos)
        annotate_text.append(cat_vals)

    field_labels = sep.join(sorted_cat)
    field_colors = sep.join([cat2style.get(c) for c in sorted_cat])
    if legend:
        l_text = deduced_legend(cat2style, sep=sep, info_name=dataset_name)
    else:
        l_text = '' 
    template_text = template_text.format(
        dataset_label=dataset_name,
        field_colors=field_colors,
        field_labels=field_labels,
        legend_text=l_text,
    )
    final_text = template_text + "\n" + "\n".join(annotate_text)
    return final_text
