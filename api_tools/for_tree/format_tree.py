from os.path import join,dirname,exists

from ete3 import Tree

indir = '/home-user/thliao/template_txt/'
if not exists(indir):
    indir = join(dirname(__file__), 'itol_template')
    
dataset_symbol_template = join(indir, 'dataset_symbols_template.txt')

def read_tree(in_tree, format=None):
    if isinstance(in_tree, str):
        t = Tree(open(in_tree).read(), format=format)
    elif isinstance(in_tree, Tree):
        t = in_tree
    else:
        raise IOError('unknown input')
    return t


def earse_name(in_tree_file, format=0):
    t = read_tree(in_tree_file, format=format)
    for n in t.traverse():
        if not n.is_leaf():
            n.name = ''
    return t


def root_tree_with(in_tree_file, gene_names=[], format=0):
    leaf_list = []
    t = read_tree(in_tree_file, format=format)
    all_leafs = t.get_leaves()
    if not gene_names:
        return t

    for gname in gene_names:
        leafs = [_ for _ in all_leafs if gname in _.name]
        leaf_list += leafs
    if len(leaf_list) != 1:
        LCA = t.get_common_ancestor(leaf_list)
        t.set_outgroup(LCA)
    elif len(leaf_list) == 1:
        t.set_outgroup(leaf_list[0])
    else:
        print("No leaf could found with input '%s'" % str(gene_names))
        return t
    return t


def sort_tree(in_tree_file, ascending=True, format=0):
    # from bottom to top
    # sort_by_num of nodes
    # ascending is True, mean branch have less leafs place bottom.
    t = read_tree(in_tree_file, format=format)
    for n in t.traverse():
        # childrens = n.children
        # if len(childrens)==2:
        # d1,d2 = [len(_.get_leaves()) for _ in n.children]
        sort_by_ascending = list(sorted(n.children,
                                        key=lambda x: len(x.get_leaves()),
                                        ))
        if ascending:
            n.children = sort_by_ascending[::-1]
        else:
            n.children = sort_by_ascending
    return t


def renamed_tree(in_tree_file, outfile=None, format=0):
    count = 0
    t = read_tree(in_tree_file, format=format)
    # t = sort_tree(in_tree_file,ascending=ascending,format=format)
    for n in t.traverse():
        if not n.name:
            if not str(int(n.support)) or str(int(n.support)) == '1.0':
                n.name = 'I%s' % (count)
            else:
                n.name = 'I%s_S%s' % (count, str(int(n.support)))
            count += 1
        elif isinstance(n.name,str) and n.name.startswith('I'):
            S_ori = n.name.split('_')[-1]
            S_ori = S_ori.strip('S')
            if not S_ori.isnumeric():
                #print(S_ori)
                S_ori = '100'
            n.name = f'I{count}_S{S_ori}'
            count +=1
    if outfile is None:
        return t
    t.write(outfile=outfile, format=3)
    return t


def add_cal_api(in_tree_file, out_newick, calibration_txt, format=0):
    t = read_tree(in_tree_file, format=format)
    calibration_dict = {}
    for row in open(calibration_txt):
        if row and not row.startswith('#'):
            LCA, time, _remained = row.split('\t')
            calibration_dict[LCA] = time
    t = earse_name(t)
    for LCA, time in calibration_dict.items():
        names = LCA.split('|')
        LCA_node = t.get_common_ancestor([l
                                          for l in t.get_leaves()
                                          if l.name in names])
        LCA_node.name = time
    final_tree = t.write(format=8, format_root_node=True).replace('NoName', '')

    text = f'{len(t.get_leaves())}\n' + final_tree
    with open(out_newick, 'w') as f1:
        f1.write(text)
    return text

def draw_cal_itol(calibration_txt,odir):
    itol_odir = odir
    size = '12'
    shape = '2'
    filled = '1'
    color = '#ff9900'

    rows_str = []
    rows = []
    for row in open(calibration_txt):
        if row and not row.startswith('#'):
            LCA, time, _remained = row.split('\t')

            row = '\t'.join([LCA, shape, size, color, filled, '1', time])
            rows.append(row)
            row = '\t'.join([LCA, time, '1', '#FF0000', 'bold', '2', '0'])
            rows_str.append(row)
    template_text = open(dataset_symbol_template).read()
    annotate_text = '\n'.join(rows)
    template_text = template_text.format(dataset_label='calibration',
                                            legend_text='',
                                            maximum_size=size)
    with open(join(itol_odir, 'dating_tree_calibration.txt'), 'w') as f1:
        f1.write(template_text + '\n' + annotate_text)

    template_text = open('/home-user/thliao/template_txt/dataset_text_template.txt').read()
    annotate_text = '\n'.join(rows_str)
    with open(join(itol_odir, 'dating_tree_calibration_str.txt'), 'w') as f1:
        f1.write(template_text + '\n' + annotate_text)
