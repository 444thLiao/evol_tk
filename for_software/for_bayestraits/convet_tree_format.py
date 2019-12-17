import sys

from ete3 import Tree

intree = './trees/iqtree/over20p_bac120.formatted.newick'
# intree = './trees/iqtree/over20p_bac120.ufboot'
otree = './bayesTraits_test/test.trees'

intree, otree = sys.argv[1:]

root_with = 'GCA_900097105.1,GCA_000020225.1,GCA_000172155.1,GCA_001318295.1,GCA_001613545.1,GCA_000019665.1,GCA_000019965.1,GCA_001746835.1'
if __name__ == "__main__":
    if len(open(intree).read().split('\n')) == 1:
        t = Tree(intree, format=3)
    else:
        multiple_trees = []
        for row in open(intree):
            row = row.strip('\n')
            multiple_trees.append(Tree(row))

        LCA = t.get_common_ancestor(root_with.split(','))
        t.set_outgroup(LCA)
        # TODO: finish it.

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
    with open(otree, 'w') as f1:
        f1.write(final_text)
