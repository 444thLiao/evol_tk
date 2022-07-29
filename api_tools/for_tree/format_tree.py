from os.path import join, dirname, exists

from ete3 import Tree

indir = '/home-user/thliao/template_txt/'
if not exists(indir):
    indir = join(dirname(__file__), 'itol_template')

dataset_symbol_template = join(indir, 'dataset_symbols_template.txt')
dataset_text_template = join(indir, 'dataset_text_template.txt')


def read_tree(in_tree, format=None):
    if isinstance(in_tree, str) and exists(in_tree):
        if format=='auto':
            for f in [0,1,2,3,4,5]:
                try:
                    t = Tree(in_tree, format=f)
                    return t
                except:
                    pass
        else:
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
        if LCA is t:
            print("same root as previous.")
            pass
        else:
            t.set_outgroup(LCA)
    elif len(leaf_list) == 1:
        t.set_outgroup(leaf_list[0])
    else:
        print("No leaf could found with input '%s'" % str(gene_names))
        return t
    return t


def sort_tree(in_tree_file, ascending=True, format=0):
    # the order of the number of descendent leaves from bottom to top
    # sort each internal node by the number of leafs descendent of it
    # ascending is True, mean branch have less leafs place bottom.
    t = read_tree(in_tree_file, format=format)
    for n in t.traverse():
        sort_by_ascending = list(sorted(n.children,
                                        key=lambda x: len(x.get_leaves()), ))
        # default is ascending according to description of builtin function `sorted`

        if ascending:
            # due to from bottom to top instead of reverse... so it need to reverse by [::-1]
            n.children = sort_by_ascending[::-1]
        else:
            n.children = sort_by_ascending
    return t


def renamed_tree(in_tree_file, outfile=None, format=0):
    """
    rename internal nodes by bootstrap values and the index of it.

    :param in_tree_file:
    :param outfile:
    :param format:
    :return:
    """
    count = 0
    t = read_tree(in_tree_file, format=format)
    for n in t.traverse():
        if not n.name:
            if not str(int(n.support)) or str(int(n.support)) == '1':
                n.name = 'I%s' % (count)
            else:
                n.name = 'I%s_S%s' % (count, str(int(n.support)))
            count += 1
        elif isinstance(n.name, str) and (n.name.startswith('I') and '_' in n.name ):
            S_ori = n.name.split('_')[-1]
            S_ori = S_ori.strip('S')
            if not S_ori.isnumeric():
                # print(S_ori)
                S_ori = '100'
            n.name = f'I{count}_S{S_ori}'
            count += 1
    if outfile is None:
        return t

    return t.write(outfile=outfile, format=3)


def add_cal_api(in_tree_file, out_newick, calibration_txt, format=0):
    """
    Implement a function to add calibration into tree in MCMCTree format.

    :param in_tree_file:
    :param out_newick:
    :param calibration_txt:
    :param format:
    :return:
    """
    t = read_tree(in_tree_file, format=format)
    calibration_dict = {}
    # iterate all rows of input calibration txt
    # stodge the information into calibration_dict
    for row in open(calibration_txt):
        if row.strip('\n') and not row.startswith('#'):
            LCA, time, _remained = row.split('\t')
            calibration_dict[LCA] = time
    # clean all internal names of input tree
    t = earse_name(t)
    # iterate each calibration in calibration_dict
    for LCA, time in calibration_dict.items():
        if LCA.upper() in ['ROOT','OROOT'] :
            LCA_node = t
        else:
            names = LCA.split('|')
            # get the common ancestor
            LCA_node = t.get_common_ancestor([l
                                            for l in t.get_leaves()
                                            if l.name in names])
        # rename the common ancestor with give name
        LCA_node.name = time
    # write out the tree
    final_tree = t.write(format=8, format_root_node=True).replace('NoName', '')
    # format the time information into suitable format
    for v in calibration_dict.values():
        if '(' in v:
            final_tree = final_tree.replace(v.replace('(', '_').replace(')', '_').replace(',', '_'),
                                            "'%s'" % v)
    text = f'{len(t.get_leaves())}\n' + final_tree
    with open(out_newick, 'w') as f1:
        f1.write(text)
    return text


def draw_cal_itol(calibration_txt, odir):
    """
    abbr: draw calibration nodes and text in itol required format

    Output files is force named 'dating_tree_calibration.txt' and 'dating_tree_calibration_str.txt'
    :param calibration_txt:
    :param odir:
    :return:
    """
    itol_odir = odir
    size = '12'
    shape = '2'
    filled = '1'
    color = '#ff9900'

    rows_str = []
    rows = []
    for row in open(calibration_txt):
        if row and not row.startswith('#'):
            try:
                LCA, time, _remained = row.split('\t')
            except:
                LCA, time = row.split('\t')

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

    template_text = open(dataset_text_template).read()
    annotate_text = '\n'.join(rows_str)
    with open(join(itol_odir, 'dating_tree_calibration_str.txt'), 'w') as f1:
        f1.write(template_text + '\n' + annotate_text)
