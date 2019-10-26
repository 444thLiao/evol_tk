from ete3 import Tree
def root_tree_with(in_tree_file,gene_names=[],format=0):
    leaf_list = []
    t = Tree(open(in_tree_file).read(),format=format)
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
    if isinstance(in_tree_file,Tree):
        t = in_tree_file
    else:
        t = Tree(open(in_tree_file).read(),format=format)
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