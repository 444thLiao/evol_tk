from collections import defaultdict
from os.path import *


def _get_tophit(gid2locus, top_hit):
    if top_hit:
        gid2locus = {k: sorted(v,
                               key=lambda x: x[1])
                     for k, v in gid2locus.items()}
        gid2locus = {k: [v[0][0]]
        if v else []
                     for k, v in gid2locus.items()}
    else:
        gid2locus = {k: [_[0] for _ in v]
        if v else []
                     for k, v in gid2locus.items()}
    return gid2locus


def _parse_blastp(ofile, match_ids=[], filter_evalue=1e-3, top_hit=False):
    if not match_ids:
        gid2locus = defaultdict(list)
    else:
        gid2locus = {k: [] for k in match_ids}
    for row in open(ofile, 'r'):
        sep_v = row.split('\t')
        locus = sep_v[0]
        evalue = float(sep_v[10])
        if filter_evalue and evalue > filter_evalue:
            continue
        if sep_v[1] in match_ids:
            gid2locus[sep_v[1]].append((locus, evalue))
        if not match_ids:
            gid2locus[sep_v[1]].append((locus, evalue))
    gid2locus = _get_tophit(gid2locus, top_hit=top_hit)
    return gid2locus


def _parse_hmmscan(ofile, filter_evalue=1e-20, top_hit=False, gene_pos=0):
    gid2locus = defaultdict(list)

    for row in open(ofile, 'r'):
        if row.startswith('#'):
            continue
        r = row.split(' ')
        r = [_ for _ in r if _]

        gene_id = r[gene_pos]
        locus_tag = r[2]
        evalue = float(r[4])
        if filter_evalue and evalue > filter_evalue:
            continue
        gid2locus[gene_id].append((locus_tag, evalue))
    if top_hit:
        gid2locus = _get_tophit(gid2locus, top_hit=top_hit)
    return gid2locus


def process_path(path):
    if '~' in path:
        path = expanduser('path')
    if not '/' in path:
        path = './' + path
    path = abspath(path)
    return path


# two function for dating workflow(formatting the name of assembly genomes)
def convert_genome_ID(genome_ID):
    # for GCA_900078535.2
    # it will return
    if isinstance(genome_ID, str) and genome_ID.startswith('GC'):
        return genome_ID.split('_')[-1].replace('.', 'v')
    else:
        return genome_ID


def convert_genome_ID_rev(locus_ID):
    # for 900078535v2
    # it will return
    if isinstance(locus_ID, str):
        if '|' in locus_ID:
            # other labmater used
            genome_name = locus_ID.partition('|')[0]
            return genome_name
        if '_' in locus_ID:
            # tianhua version, it won't contain |
            locus_ID = locus_ID.partition('_')[0]
            return 'GCA_' + locus_ID.replace('v', '.')
        else:
            return 'GCA_' + locus_ID.replace('v', '.')
    else:
        return locus_ID
