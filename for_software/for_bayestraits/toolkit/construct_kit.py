from ete3 import Tree


# intree = '../trees/iqtree/over20p_bac120.formatted.newick'
# intree = './trees/iqtree/over20p_bac120.ufboot'
# ocommand = './tag.command'

def get_tags(intree):
    """
    generate cmd for 'ADD all internal node'
    Use formatted newick which generate by thliao. (internal node has been named one)
    """
    rows = []
    t = Tree(intree, format=1)
    for _ in t.traverse():
        if (not _.is_leaf()) and (_.name):
            rows.append(f"AddTag {_.name} {' '.join(_.get_leaf_names())}")
            rows.append(f"AddNode {_.name} {_.name}")
    return rows


def nw2nexus(t, root_with=''):
    new_name2old_name = {}
    _count = 0
    for leaf in t.get_leaves():
        new_name2old_name[str(_count)] = leaf.name
        leaf.name = str(_count)
        _count += 1

    # for bayestraits
    nexus_template = """#NEXUS
    begin trees;
                translate
{translate_text};

                        tree tree1 = {tree_text}
                        
    end;
    """

    translate_text = ',\n'.join([20 * ' ' + f"{c} {v}" for c, v in new_name2old_name.items()])
    tree_text = t.write(format=5)

    final_text = nexus_template.format(translate_text=translate_text,
                                       tree_text=tree_text)
    return final_text
