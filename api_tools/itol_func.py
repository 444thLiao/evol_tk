from ete3 import Tree
from os.path import join,exists,dirname
import plotly.express as px

indir = '/home-user/thliao/template_txt/'
if not exists(indir):
    indir = join(dirname(__file__),'itol_template')

color_strip_template = join(indir, 'dataset_color_strip_template.txt')
dataset_styles_template = join(indir, 'dataset_styles_template.txt')
dataset_binary_template = join(indir, 'dataset_binary_template.txt')
label_template = join(indir,'labels_template.txt')
dataset_symbol_template = join(indir,'dataset_symbols_template.txt')
matrix_like_template = join(indir,"dataset_external_shapes_template.txt")


def root_tree_with(in_tree_file,gene_names=[],format=0):
    leaf_list = []
    t = Tree(open(in_tree_file).read(),format=format)
    all_leafs = t.get_leaves()
    if not gene_names:
        return t
    for gname in gene_names:
        leafs = [_ for _ in all_leafs if gname in _.name]
        leaf_list+= leafs
    LCA = t.get_common_ancestor(leaf_list)
    t.set_outgroup(LCA)
    return t

def sort_tree(in_tree_file,ascending=True,format=0):
    # from top to bottom
    # ascending is True, mean longer branch place bottom.
    if isinstance(in_tree_file,Tree):
        t = in_tree_file
    else:
        t = Tree(open(in_tree_file).read(),format=format)
    for n in t.traverse():
        childrens = n.children
        if len(childrens)==2:
            d1,d2 = [_.dist for _ in n.children]
            if ascending:
                if d1>d2:
                    n.children = n.children[::-1]
                elif d1==d2:
                    # place outgroup at the bottom
                    n.children = list(sorted(n.children,
                                         key=lambda x:len(x.get_leaves()),
                                         ))[::-1]
            else:
                if d1<d2:
                    n.children = n.children[::-1]
                elif d1==d2:
                    # place outgroup at the top
                    n.children = list(sorted(n.children,
                                         key=lambda x:len(x.get_leaves()),
                                         ))
    return t

def renamed_tree(in_tree_file, outfile,ascending=True,format=0):
    count = 0
    #t = Tree(open(in_tree_file).read())
    t = sort_tree(in_tree_file,ascending=ascending,format=format)
    for n in t.traverse():
        if not n.name:
            n.name = 'I%s_S%s' % (count,str(int(n.support)))
            count += 1
    t.write(outfile=outfile, format=3)
    return t


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


def deduced_field(info2style, infos, sep='\t'):
    template_text = open(dataset_binary_template).read()
    colors_theme = px.colors.qualitative.Light24
    shapes = []
    labels = []
    colors = []
    for idx, info in enumerate(infos):
        shapes.append(info2style[info].get('shape', '3'))
        labels.append(info2style[info].get('info', info))
        colors.append(info2style[info].get('color', colors_theme[idx]))
    template_text = template_text.format(field_shapes=sep.join(shapes),
                                         field_labels=sep.join(labels),
                                         field_colors=sep.join(colors))
    return template_text


def annotate_outgroup(ID2infos, info2style,):
    for ID, infos in ID2infos.items():
        # arbitary to choose a column of infos as template for deduced legend
        break
    template_text = deduced_field(info2style, infos)
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

def to_color_strip(ID2info, info2color, info_name='dataset'):
    template_text = open(color_strip_template).read()
    id2col = {id: info2color[info] for id, info in ID2info.items()}
    annotate_text = '\n'.join(['%s,%s\n' % (id, col)
                               for id, col in id2col.items()])
    legend_text = deduced_legend(info2color, info_name)

    template_text = template_text.format(legend_text=legend_text,
                                         dataset_label=info_name)
    info_name = info_name.replace('/', '_')
    return template_text+'\n'+annotate_text


def to_color_labels_bg(ID2info, info2color, info_name='labels bg'):
    # clade for
    template_text = open(dataset_styles_template).read()
    id2col = {ID: info2color[info] for ID, info in ID2info.items()}
    each_template = '{ID}\t{TYPE}\t{WHAT}\t{COLOR}\t{WIDTH_OR_SIZE_FACTOR}\t{STYLE}\t{BACKGROUND_COLOR}\n'
    #legend_text = deduced_legend(info2color, info_name, sep='\t')

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
    #normally for bootstrap
    # give it a tree is enough, must have internal node name.
    template_text = open(dataset_symbol_template).read()
    tree = Tree(in_tree,format=1)
    id2support = {}
    for n in tree.traverse():
        if (not n.is_leaf()) and n.name:
            support_v = int(n.name.split('_S')[-1])
            id2support[n.name] = support_v
    
    # ID,symbol,size,color,fill,position,label
    rows = []
    for id,s_v in id2support.items():
        size = '5'
        shape = '2'
        filled = '1'
        if int(s_v) >= 95:
            color = '#000000'
        else:
            color = '#999999'
        row = '\t'.join([id,shape,size,color,filled,'0',''])
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
    template_text = template_text.format(dataset_label=dataset_label,
                                         field_color=color,
                                         field_labels='\t'.join(all_v))
    annotate_text = ''
    for ID,v in ID2categorical_v.items():
        vals = [ID]
        for _ in range(len(all_v)):
            if v == all_v[_]:
                vals.append('10')
            else:
                vals.append('0')
        annotate_text += '\t'.join(vals) + '\n'
    return template_text + annotate_text