from Bio.Data import IUPACData
from Bio.Seq import Seq
import re
tmp_d = IUPACData.ambiguous_dna_values
def get_region(seq,f,r):
    # use a pair of primer to retrieve the amplified region
    
    f = ''.join([l if len(tmp_d.get(l,l))==1 else f"[{tmp_d.get(l,l)}]" for l in f])
    r = str(Seq(r).reverse_complement())
    f = ''.join([l if len(tmp_d.get(l,l))==1 else f"[{tmp_d.get(l,l)}]" for l in f])
    r = ''.join([l if len(tmp_d.get(l,l))==1 else f"[{tmp_d.get(l,l)}]" for l in r])
    region_include_primers = re.search(f+'([ATCG]*)'+r,seq)
    if region_include_primers is None:
        return 'No match'
    # return the sequence removing the primer
    return region_include_primers.groups()[0]
    