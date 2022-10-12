
import random
from collections import defaultdict
def denovo_cluster(tre,leaves,sampled=True):
    name2node = {_.name:_ for _ in tre.traverse()}
    assert len(name2node) == len(list(tre.traverse())) 
    # make sure no nodes with identical name
    final_nodes = leaves[::]
    _c = 0
    while 1:
        #print(f"{_c} times")
        ori_len = len(final_nodes)
        p2child = defaultdict(list)
        for n in final_nodes:
            p2child[name2node[n].up.name].append(n)
        new_final = []
        for k,v in p2child.items():
            if len(v)>1:
                new_final.append(k)
            elif len(v)==1:
                new_final.append(v[0])
            else:
                raise IOError()
        new_len = len(new_final)
        if ori_len == new_len:
            break
        final_nodes = new_final[::]
        _c+=1
    if not sampled:
        return final_nodes
    else:
        sampled_gids = []
        for n in final_nodes:
            if name2node[n].is_leaf():
                sampled_gids.append(n)
            else:
                sampled_gids.append(random.choice(name2node[n].get_leaf_names()))
        return sampled_gids