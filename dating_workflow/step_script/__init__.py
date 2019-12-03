
from collections import defaultdict
def _get_tophit(gid2locus,top_hit):
                
    if top_hit:
        gid2locus = {k:sorted(v,
                            key=lambda x:x[1])  
                      for k,v in gid2locus.items()}
        gid2locus = {k:[v[0][0]] 
                      if v else [] 
                      for k,v in gid2locus.items()}
    else:
        gid2locus = {k:[_[0] for _ in v] 
                      if v else []
                      for k,v in gid2locus.items()}
    return gid2locus

def _parse_blastp(ofile,match_ids=[],top_hit = False):
    if not match_ids:
        gid2locus = defaultdict(list)
    else:
        gid2locus = {k:[] for k in match_ids}
    for row in open(ofile,'r'):
        sep_v = row.split('\t')
        locus = sep_v[0]
        evalue = sep_v[10]
        if sep_v[1] in match_ids:
            gid2locus[sep_v[1]].append((locus,evalue))
        if not match_ids:
            gid2locus[sep_v[1]].append((locus,evalue))
    gid2locus = _get_tophit(gid2locus,top_hit=top_hit)
    return gid2locus

def _parse_hmmscan(ofile,filter_evalue=None,top_hit = False):
    gid2locus = defaultdict(list)

    for row in open(ofile, 'r'):
        if row.startswith('#'):
            continue
        r = row.split(' ')
        r = [_ for _ in r if _]

        gene_id = r[1]
        locus_tag = r[2]
        evalue = float(r[4])
        if filter_evalue and evalue <= filter_evalue:
            gid2locus[gene_id].append((locus_tag, evalue))
        else:
            gid2locus[gene_id].append((locus_tag, evalue))
    gid2locus = _get_tophit(gid2locus,top_hit=top_hit)
    return gid2locus
