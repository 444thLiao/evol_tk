"""
Automatic sampling according defined criterias for following dating analysis

It would consider
1. the completeness of target genes (now is cog25)
2. the number of descending leaves to the total number of genomes it sampled
3. is it a basal one
"""


from collections import defaultdict
from email.policy import default
from pickle import NONE


def precluster_based_selection(nodes, l2cluster=None):
    if l2cluster is None:
        return nodes
    haved_c = []
    remained_nodes = []
    for n in nodes:
        if l2cluster[n] not in haved_c:
            remained_nodes.append(n)
            haved_c.append(l2cluster[n])
    return remained_nodes


def get_cluster(c_file):
    cluster2genomes = defaultdict(list)
    for row in open(c_file).read().split('\n')[1:]:
        if not row:
            continue
        rows = row.split('\t')
        if rows[1] != '-1':
            cluster2genomes[rows[1]].append(rows[0])
        else:
            cluster2genomes[rows[0]].append(rows[0])
    return cluster2genomes


def get_basal_ones(node, up_level=3, maximum_retained=3, **kwargs):
    # get required genoems which could maintain the topological structures of request node up to the root.
    all_leaves = []
    curr_n = node
    cur_l = 0
    while cur_l <= up_level:
        if curr_n.is_root():
            break
        sister_node = [_
                       for _ in curr_n.up.children
                       if _.name != curr_n.name][0]
        leaves = get_simple_LCA(
            sister_node, maximum=maximum_retained, **kwargs)
        all_leaves += leaves
        curr_n = curr_n.up
        cur_l += 1
        # print(len(all_leaves))
    return all_leaves


def get_simple_LCA(node, maximum=10, l2cluster=None,
                   validated_func=lambda x: True,
                   sort_func=None
                   ):
    # get the necessary genomes descending from the requested node
    match_leaves = [n
                    for n in node.get_leaf_names()
                    if validated_func(n)]

    l2dis = {n.name: node.get_distance(n, topology_only=True)
             for n in node.get_leaves()}
    l2dis = {k: v for k, v in l2dis.items() if k in match_leaves}
    
    if sort_func is None:
        def sort_func(x): return sorted(x, key=lambda g: l2dis[g])
    if len(match_leaves) <= maximum:
        if len(match_leaves)==0:
            raise IOError
        return match_leaves
    else:
        sorted_leaves = sort_func(match_leaves)
        sorted_leaves = precluster_based_selection(
            sorted_leaves, l2cluster)[:maximum]
        if len(sorted_leaves)==0:
            raise IOError
        return sorted_leaves


def sample_children(tree_node):
    while 1:
        smaller_n = sorted(tree_node.children,
                           key=lambda x: len(x.get_leaf_names()))[0]


def sampling(st, target_nodes_text,
             must_in=[], node2cluster=None,
             up_level=2, max_num_up_each_level=3,
             sort_func=None,
             max_num_down=10, validated_func=lambda x: True):
    """
    node2cluster: leaf 2 cluster name

    sample_l2reasons has NN kinds of results:
    given: passed to the parameter 'must_in'
    Descendent
    Outgroup
    """
    sample_l2reasons = defaultdict(list)
    [sample_l2reasons[_].append('given') for _ in must_in]
    st = st.copy()
    name2nodes = {n.name: n for n in st.traverse()}
    target_nodes = [name2nodes[_] for _ in target_nodes_text]

    for tn in target_nodes:
        leaves1 = []
        for child in tn.children:
            leaves1 += get_simple_LCA(child,
                                      maximum=max_num_down//2,
                                      l2cluster=node2cluster, validated_func=validated_func,
                                      sort_func=sort_func)

        [sample_l2reasons[_].append(f'Descendant of {tn.name}')
         for _ in leaves1]
        leaves2 = get_basal_ones(tn,
                                 up_level=up_level,
                                 maximum_retained=max_num_up_each_level,
                                 l2cluster=node2cluster,
                                 validated_func=validated_func,
                                 sort_func=sort_func)
        [sample_l2reasons[_].append(f'Outgroup of {tn.name}')
         for _ in leaves2]

    return sample_l2reasons


if __name__ == "__main__":
    pass
