import sys
from os.path import join,exists,dirname
sys.path.insert(0,dirname(__file__))
from ete3 import Tree
import plotly.express as px
from .for_tree.format_tree import *

indir = '/home-user/thliao/template_txt/'
if not exists(indir):
    indir = join(dirname(__file__),'itol_template')

color_strip_template = join(indir, 'dataset_color_strip_template.txt')
dataset_styles_template = join(indir, 'dataset_styles_template.txt')
dataset_binary_template = join(indir, 'dataset_binary_template.txt')
label_template = join(indir,'labels_template.txt')
dataset_symbol_template = join(indir,'dataset_symbols_template.txt')
matrix_like_template = join(indir,"dataset_external_shapes_template.txt")
labels_template = join(indir,"labels_template.txt")
dataset_text_template = join(indir,"dataset_text_template.txt")
def deduced_legend(info2color, info_name='dataset', sep=','):
    # for implemented a legend with dictinonary named info2color.
    
    legend_title = info_name
    legend_shape = sep.join(['1'] * len(info2color))
    legend_colors = sep.join([_
                              for _ in info2color.values()])
    legend_labels = sep.join(list(info2color.keys()))
    legend_text = f"""
LEGEND_TITLE{sep}{legend_title}
LEGEND_SHAPES{sep}{legend_shape}
LEGEND_COLORS{sep}{legend_colors}
LEGEND_LABELS{sep}{legend_labels}"""
    return legend_text


def deduced_legend2(info2style, infos, sep='\t'):
    # for info2style instead of info2color
    colors_theme = px.colors.qualitative.Dark24 + px.colors.qualitative.Light24
    shapes = []
    labels = []
    colors = []
    for idx, info in enumerate(infos):
        shapes.append(info2style[info].get('shape', '1'))
        labels.append(info2style[info].get('info', info))
        colors.append(info2style[info].get('color', colors_theme[idx]))
    legend_text = ['FIELD_SHAPES'+sep + sep.join(shapes),
                   'FIELD_LABELS'+sep + sep.join(labels),
                   'FIELD_COLORS'+sep + sep.join(colors),]
    legend_text = '\n'.join(legend_text)
    return legend_text

def annotate_outgroup(ID2infos, info2style,):
    for ID, infos in ID2infos.items():
        # arbitary to choose a column of infos as template for deduced legend
        break
    template_text = deduced_legend2(info2style, infos)
    annotate_text = ''
    for ID, infos in ID2infos.items():
        row = [ID]
        for info in infos:
            if info not in info2style:
                row.append('-1')
            elif info in info2style:
                row.append(info2style[info]['status'])
        annotate_text += '\t'.join(row) + '\n'
    return template_text + annotate_text

def to_binary_shape(ID2info,info2style=None, info_name='dataset',manual_v=[],omitted_other=False):
    # id2info, could be {ID:list/set}
    # info2color: could be {gene1: {shape:square,color:blabla},}
    # None will use default.

    template_text = open(dataset_binary_template).read()
    if not manual_v:
        all_v = list(sorted(set([_ for v in ID2info.values() for _ in v if _])))
    else:
        all_v = manual_v
    if info2style is None:
        info2style = {k:{} for k in all_v}
    others_label = '-1' if omitted_other else '0'
    annotate_text = []
    for ID,vset in ID2info.items():
        row = '\t'.join([ID] + ['1' if _ in vset else others_label for _ in all_v])
        annotate_text.append(row)
    annotate_text = '\n'.join(annotate_text)
    
    legend_text = deduced_legend2(info2style,all_v,sep='\t')
    template_text = template_text.format(legend_text=legend_text,
                                         dataset_label=info_name)
    return template_text+'\n'+annotate_text

def to_color_strip(ID2info, info2color, info_name='dataset'):
    template_text = open(color_strip_template).read()
    id2col = {id: info2color[info] for id, info in ID2info.items()}
    annotate_text = '\n'.join(['%s,%s\n' % (id, col)
                               for id, col in id2col.items()])
    legend_text = deduced_legend(info2color, info_name)
    info_name = info_name.replace('/', '_')
    template_text = template_text.format(legend_text=legend_text,
                                         dataset_label=info_name)
    
    return template_text+'\n'+annotate_text


def to_color_labels_bg(ID2info, info2color, info_name='labels bg'):
    # clade for
    template_text = open(dataset_styles_template).read()
    id2col = {ID: info2color[info] for ID, info in ID2info.items()}
    each_template = '{ID}\t{TYPE}\t{WHAT}\t{COLOR}\t{WIDTH_OR_SIZE_FACTOR}\t{STYLE}\t{BACKGROUND_COLOR}\n'
    legend_text = deduced_legend(info2color, info_name, sep='\t')

    template_text = template_text.format(dataset_label=info_name,
                                         legend_text='')
    
    rows = [each_template.format(ID=ID,
                                  TYPE='label',
                                  WHAT='node',
                                  COLOR='#000000',
                                  WIDTH_OR_SIZE_FACTOR=1,
                                  STYLE='bold',
                                  BACKGROUND_COLOR=color)
             for ID, color in id2col.items()]
    return template_text + '\n'.join(rows)

def to_color_branch(ID2info, info2color, dataset_name='color branch',no_legend=False):
    # clade for
    template_text = open(dataset_styles_template).read()
    id2col = {ID: info2color[info] for ID, info in ID2info.items()}
    each_template = '{ID}\t{TYPE}\t{WHAT}\t{COLOR}\t{WIDTH_OR_SIZE_FACTOR}\t{STYLE}\t{BACKGROUND_COLOR}\n'
    if no_legend:
        legend_text = ''
    else:
        legend_text = deduced_legend(info2color, dataset_name, sep='\t')

    template_text = template_text.format(dataset_label=dataset_name,
                                         legend_text=legend_text)
    rows = [each_template.format(ID=ID,
                                 TYPE='branch',
                                 WHAT='node',
                                 COLOR=color,
                                 WIDTH_OR_SIZE_FACTOR=3,
                                 STYLE='normal',
                                 BACKGROUND_COLOR='')
            for ID, color in id2col.items()]
    rows += [each_template.format(ID=ID,
                                  TYPE='label',
                                  WHAT='node',
                                  COLOR=color,
                                  WIDTH_OR_SIZE_FACTOR=1,
                                  STYLE='bold',
                                  BACKGROUND_COLOR='')
             for ID, color in id2col.items()]
    return template_text + '\n'.join(rows)


def to_color_Clade(ID2info, info2color, tree, 
                   dataset_name='color branch clade',
                   no_legend=False,bgcolor={}):
    # colorize branch within whole clade when all children under this clade share same info.
    template_text = open(dataset_styles_template).read()
    id2col = {ID: info2color[info] for ID, info in ID2info.items()}
    each_template = '{ID}\t{TYPE}\t{WHAT}\t{COLOR}\t{WIDTH_OR_SIZE_FACTOR}\t{STYLE}\t{BACKGROUND_COLOR}'
    if no_legend:
        legend_text = ''
    else:
        legend_text = deduced_legend(info2color, dataset_name, sep='\t')

    tree_obj = Tree(tree, format=3)

    def collapsed_leaf(node):
        if node.is_leaf():
            return True
        else:
            leafs_all = node.get_leaves()
            children_names = [ID2info.get(_.name, '') for _ in leafs_all]
            if '' in children_names:
                return False
            if len(set(children_names)) == 1:
                return True
            else:
                return False
        return False
    new_tree_obj = Tree(tree_obj.write(is_leaf_fn=collapsed_leaf))
    new_leaves_names = [_.name for _ in new_tree_obj.get_leaves()]
    internal_nodes = [tree_obj.search_nodes(
        name=_)[0] for _ in new_leaves_names]
    internal_nodes = [_ for _ in internal_nodes if not _.is_leaf()]
    internal_node2info = {n.name: info2color.get(ID2info.get(n.get_leaves()[0].name, ''), '')
                          for n in internal_nodes}
    internal_node2info = {n.name: info2color[ID2info[n.get_leaves()[0].name]]
                          for n in internal_nodes}
    #dropped_IDs = [_.name for node in internal_nodes for _ in node.get_leaves()]
    # for ID in list(id2col.keys()):
    #     id2col.pop(ID)

    template_text = template_text.format(dataset_label=dataset_name,
                                         legend_text=legend_text)
    rows = [each_template.format(ID=ID,
                                 TYPE='branch',
                                 WHAT='node',
                                 COLOR=color,
                                 WIDTH_OR_SIZE_FACTOR=7,
                                 STYLE='normal',
                                 BACKGROUND_COLOR='')
            for ID, color in id2col.items()]
    
    rows += [each_template.format(ID=ID,
                                  TYPE='label',
                                  WHAT='node',
                                  COLOR=color,
                                  WIDTH_OR_SIZE_FACTOR=1,
                                  STYLE='bold',
                                  BACKGROUND_COLOR=bgcolor.get(ID,''))
             for ID, color in id2col.items()]
    rows += [each_template.format(ID=ID,
                                  TYPE='branch',
                                  WHAT='clade',
                                  COLOR=color,
                                  WIDTH_OR_SIZE_FACTOR=7,
                                  STYLE='normal',
                                  BACKGROUND_COLOR='')
             for ID, color in internal_node2info.items()]
    return template_text + '\n'.join(rows)


def to_node_symbol(in_tree,dataset_name='bootstrap'):
    # normally for bootstrap
    # give it a tree is enough, must have internal node name.
    template_text = open(dataset_symbol_template).read()
    tree = Tree(in_tree,format=3)
    id2support = {}
    for n in tree.traverse():
        if (not n.is_leaf()) and n.name:
            support_v = int(n.name.split('_S')[-1])
            id2support[n.name] = support_v
    
    # ID,symbol,size,color,fill,position,label
    rows = []
    for id,s_v in id2support.items():
        size = '7'
        shape = '2'
        filled = '1'
        if int(s_v) >= 95:
            color = '#000000'
        elif int(s_v) >= 85 and int(s_v) < 95:
            color = '#777777'
        elif int(s_v) >= 65 and int(s_v) < 85:
            color = '#eeeeee'
        else:
            color = ''
        if color:
            row = '\t'.join([id,shape,size,color,filled,'1',''])
            rows.append(row)
        
    annotate_text = '\n'.join(rows)
    template_text = template_text.format(dataset_label=dataset_name,
                                         legend_text='',
                                         maximum_size=size)
    return template_text + annotate_text



def to_matrix_shape(ID2categorical_v,dataset_label,color='#000000'):
    template_text = open(matrix_like_template).read()
    all_v = set(map(str,ID2categorical_v.values()))
    all_v = list(sorted(all_v))
    if type(color) == str:
        color_str = '\t'.join([color]*len(all_v))
    elif type(color) == dict:
        color_str = '\t'.join([color.get(_,'#000000') for _ in all_v])
    else:
        raise IOError
    template_text = template_text.format(dataset_label=dataset_label,
                                         field_color=color_str,
                                         field_labels='\t'.join(all_v))
    annotate_text = ''
    for ID,v in ID2categorical_v.items():
        vals = [ID]
        for _ in range(len(all_v)):
            if v == str(all_v[_]):
                vals.append('10')
            else:
                vals.append('0')
        annotate_text += '\t'.join(vals) + '\n'
    return template_text + annotate_text

def to_label(id2new_id):
    template_text = open(labels_template).read()
    full_text = template_text[::]
    for id,new_id in id2new_id.items():
        id = str(id)
        new_id = str(new_id)
        full_text += '%s,%s\n' % (id, new_id)
    return full_text