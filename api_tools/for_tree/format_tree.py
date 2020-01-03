from ete3 import Tree


def read_tree(in_tree,format):    
    if isinstance(in_tree,str):
        t = Tree(open(in_tree).read(),format=format)
    elif isinstance(in_tree,Tree):
        t = in_tree
    else:
        raise IOError('unknown input')
    return t

def earse_name(in_tree_file):
    t = read_tree(in_tree_file)
    for n in t.traverse():
        if not n.is_leaf():
            n.name = ''
    return t
def root_tree_with(in_tree_file,gene_names=[],format=0):
    leaf_list = []
    t = read_tree(in_tree_file,format=format)
    all_leafs = t.get_leaves()
    if not gene_names:
        return t

    for gname in gene_names:
        leafs = [_ for _ in all_leafs if gname in _.name]
        leaf_list+= leafs
    if len(leaf_list) !=1:
        LCA = t.get_common_ancestor(leaf_list)
        t.set_outgroup(LCA)
    elif len(leaf_list) ==1:
        t.set_outgroup(leaf_list[0])
    else:
        print("No leaf could found with input '%s'" % str(gene_names))
        return t
    return t

def sort_tree(in_tree_file,ascending=True,format=0):
    # from bottom to top
    # sort_by_num of nodes
    # ascending is True, mean branch have less leafs place bottom.
    t = read_tree(in_tree_file,format=format)
    for n in t.traverse():
        #childrens = n.children
        #if len(childrens)==2:
            #d1,d2 = [len(_.get_leaves()) for _ in n.children]
        sort_by_ascending = list(sorted(n.children,
                                        key=lambda x:len(x.get_leaves()),
                                        ))
        if ascending:
            n.children = sort_by_ascending[::-1]
        else:
            n.children = sort_by_ascending
    return t

def renamed_tree(in_tree_file, outfile=None,ascending=True,format=0,force=False):
    count = 0
    #t = Tree(open(in_tree_file).read())
    t = sort_tree(in_tree_file,ascending=ascending,format=format)
    for n in t.traverse():
        if not n.name or force:
            n.name = 'I%s_S%s' % (count,str(int(n.support)))
            count += 1
        
    if outfile is None:
        return t
    t.write(outfile=outfile, format=3)
    return t