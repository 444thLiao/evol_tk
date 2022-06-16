from collections import defaultdict
from os.path import *
import pandas as pd


def read_table(infile,**kwargs):
    import pandas as pd
    if not exists(infile):
        exit(f"{infile} doesn't exist")
    if infile.endswith('xlsx') or infile.endswith('xls'):
        df = pd.read_excel(infile,**kwargs)
    else:
        df = pd.read_csv(infile,**kwargs)
    return df
    
def get_tophit(gid2locus,top_hit):
                
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



def parse_hmmscan_domtblout(ofile,filter_evalue=None,top_hit = False):
    # deal with --tblout
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
    gid2locus = get_tophit(gid2locus,top_hit=top_hit)
    return gid2locus

def read_summary(metadata):
    metadata_df = pd.read_csv(
        metadata, sep='\t', low_memory=False, comment='#', header=None)
    _ = open(metadata)
    header = _.readline()
    header = _.readline().strip('# \n').split('\t')
    metadata_df.columns = header
    metadata_df = metadata_df.set_index("assembly_accession")
    return metadata_df