
from ete3 import Tree
import sys
intree = '../trees/iqtree/over20p_bac120.formatted.newick'
otree = './test.trees'

intree,otree = sys.argv[1:]

if __name__ == "__main__":
    t = Tree(intree, format=3)

    new_name2old_name = {}
    _count = 0
    for leaf in t.get_leaves():
        new_name2old_name[str(_count)] = leaf.name
        leaf.name = str(_count)
        _count+=1

    # for bayestraits
    nexus_template = """#NEXUS
    begin trees;
                translate
    {translate_text};

                        tree tree1 = {tree_text}
                        
    end;
    """

    translate_text = ',\n'.join([20*' ' + f"{c} {v}" for c,v in new_name2old_name.items()])
    tree_text = t.write(format=5)

    final_text = nexus_template.format(translate_text=translate_text,
                                    tree_text=tree_text)
    with open(otree,'w') as f1:
        f1.write(final_text)

