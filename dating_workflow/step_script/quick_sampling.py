"""
Automatic sampling according defined criterias for following dating analysis

It would consider
1. the completeness of target genes (now is cog25)
2. the number of descending leaves to the total number of genomes it sampled
3. is it a basal one
"""


from collections import defaultdict


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


def check_cog25(node, genome2cog25={}):
    # if should prepare genome2cog25 first
    # from dating_workflow.step_script.extract_cog25 import parse_annotation 
    # genome2cog25 = parse_annotation(realpath('~/data/NCBI/modified_data/cog25_annotate/'), top_hit=True, evalue=1e-20)
    # the most simple way. If the most frequent one exist. Then return True
    # nonlocal genome2cog25
    if "CDD:223556" in genome2cog25.get(node, []):
        return True
    else:
        return False


def get_basal_ones(node, up_level=3, maximum_retained=3, **kwargs):
    # get required genoems which could maintain the topological structures of request node up to the root.
    all_leaves = []
    curr_n = node
    cur_l = 0
    while cur_l <= up_level:
        if curr_n.is_root():
            break
        sister_node = [_ for _ in curr_n.up.children if _.name != curr_n.name][0]
        leaves = get_simple_LCA(sister_node, maximum=maximum_retained, **kwargs)
        all_leaves += leaves
        curr_n = curr_n.up
        cur_l += 1
        # print(len(all_leaves))
    return all_leaves


def get_simple_LCA(node, maximum=10, l2cluster=None, genome2cog25={}):
    # get the necessary genomes descending from the requested node
    match_leaves = [n for n in node.get_leaf_names() if check_cog25(n, genome2cog25=genome2cog25)]
    l2dis = {n.name: node.get_distance(n, topology_only=True) for n in node.get_leaves()}
    l2dis = {k: v for k, v in l2dis.items() if k in match_leaves}
    if len(match_leaves) <= maximum:
        return match_leaves
    else:
        sorted_leaves = sorted(match_leaves, key=lambda x: l2dis[x])
        sorted_leaves = precluster_based_selection(sorted_leaves, l2cluster)[:maximum]
        return sorted_leaves


def sampling(st, target_nodes_text, 
             must_in=[], node2cluster=None,
             up_level=3, max_num_up_each_level=3,
             max_num_down=10, genome2cog25={}):
    """
    node2cluster: leaf 2 cluster name
    """
    final_leaves = list(must_in)
    st = st.copy()
    name2nodes = {n.name: n for n in st.traverse()}
    target_nodes = [name2nodes[_] for _ in target_nodes_text]

    for tn in target_nodes:
        final_leaves += get_simple_LCA(tn, maximum=max_num_down,
                                       l2cluster=node2cluster, genome2cog25=genome2cog25)
        final_leaves += get_basal_ones(tn, up_level=up_level,
                                       maximum_retained=max_num_up_each_level,
                                       l2cluster=node2cluster, genome2cog25=genome2cog25)

    final_leaves = list(set(final_leaves))
    return final_leaves


if __name__ == "__main__":
    pass
