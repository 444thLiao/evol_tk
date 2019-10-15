import pandas as pd
from os.path import *
import os
from Bio import SeqIO
from collections import defaultdict
from glob import glob
from tqdm import tqdm


blast_db = '/home-db/pub/protein_db/CDD/Cdd'
r_list = '/home-user/thliao/resource/r-prot.COG.list'
r_ids = [row.split('\t')[0] for row in open(r_list).read().split('\n') if row.split('\t')]
r_ids = [_ for _ in r_ids if _]
tbl_list = '/home-user/sswang/resource/db/cdd/cddid_all.tbl'
tbl_df = pd.read_csv(tbl_list,sep='\t',header=None,index_col=0)
# lots of id be overlapped
cdd_id2new_id = dict(zip(tbl_df.index,tbl_df.loc[:,1]))
rev_c2n = {v:k for k,v in cdd_id2new_id.items()}

cdd_id_list = [str(rev_c2n[rid]) for rid in r_ids]

def parse_tab(infile):
    sname = basename(infile).split('.')[0]
    query_df = pd.read_csv(infile,sep='\t',header=None,index_col=1)
    query_df.index = [_.split(':')[1] for _ in query_df.index]
    
    sub_cdd_id_list = [_ for _ in cdd_id_list if _ in query_df.index]
    _df = query_df.loc[sub_cdd_id_list,:]
    _df.index = [cdd_id2new_id[int(_)]
                 for _ in _df.index]
    cog2pid = dict(zip(_df.index,
                          _df.loc[:,0]))
    #pid2cog = {v:k for k,v in cog2pid.items()}
    return {sname:cog2pid}

def extract_cog(in_fas,odir,sname_dict):
    cog2records = defaultdict(list)
    for in_fa in tqdm(in_fas):
        sname = basename(in_fa).split('.')[0]
        if sname not in sname_dict:
            continue
        rid2record = {_.id:_ for _ in SeqIO.parse(in_fa,format='fasta')}
        cog2pid = sname_dict[sname]
        for cog,pid in cog2pid.items():
            record = rid2record[pid]
            record.id = sname
            cog2records[cog].append(record)
            
    if not exists(odir):
        os.makedirs(odir)
    for cog,records in cog2records.items():
        with open(join(odir,f'{cog}.faa'),'w') as f1:
            SeqIO.write(records,f1,format='fasta-2line')
            
            
if __name__ == "__main__":
    # todo: modify it into general/api version
    
    total_dict = {}
    for infile in glob('./genome_protein_files_more/cogRbp_anno/*.cogrbp'):
        total_dict.update(parse_tab(infile))
    with open('./genome_protein_files_more/extract_cog/selected_genomes.txt','w') as f1:
        required_gs = list(total_dict.keys())
        required_gs = [_ for _ in required_gs if len(total_dict[_])>=40]
        f1.write('\n'.join(required_gs))
    extract_cog(glob('./genome_protein_files_more/*.faa'),
                odir = join('./genome_protein_files_more/extract_cog'),
                sname_dict=total_dict)
    

    