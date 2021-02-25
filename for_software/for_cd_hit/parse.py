"""
parse result from cd-hit
"""
from collections import defaultdict

def parse_clstr(clstr):
    """
    make sure not special sequence id name which contains % or *
    :param clstr:
    :return:
    """
    cluster2seqs = defaultdict(list)
    cluster2repr = {}
    current_cluster = ''
    seq_name = ''
    for row in open(clstr):
        if row.startswith('>'):
            current_cluster = row.strip('>\n')
            continue
        else:
            seq_name = row.split('...')[0].split('>')[-1]
            cluster2seqs[current_cluster].append(seq_name)
        if '*' in row and not '%' in row:
            cluster2repr[current_cluster] = seq_name
    return cluster2seqs,cluster2repr

