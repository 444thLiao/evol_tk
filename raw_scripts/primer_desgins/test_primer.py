
import pandas as pd
from Bio import SeqIO
import os
from itertools import product
from Bio.Data import IUPACData
from subprocess import check_call

dgbase = {k:list(v) for k,v in IUPACData.ambiguous_dna_values.items()}

cols = ['qseqid',
 'sseqid',
 'pident',
 'length',
 'evalue',
 'mismatch',
 'gapopen',
 'sstart',
 'send',
 'stitle']

# complex version
def filtering(tab_path,fa_file,num_mismatch):
    query_len = len(list(SeqIO.parse(fa_file,'fasta'))[0].seq)
    df = pd.read_csv(tab_path,sep='\t',header=None)
    df.columns = cols
    df.loc[:,'nm'] = (query_len - df['length']*df['pident']/100)
    sub_df = df['sseqid'][df['nm']<=num_mismatch]
    return set(sub_df)

def convert_df(df):
    _tmp = df[['sseqid','sstart','send']]
    _tmp = _tmp.set_index('sseqid')
    _tmp = _tmp.loc[~_tmp.index.duplicated(),:]
    d = _tmp.to_dict('index')
    d = {k: (min([v['sstart'],v['send']]),
             max([v['sstart'],v['send']]) )  
         for k,v in d.items()}
    return d

def get_product(fp_tab,rp_tab,fp_fasta,rp_fasta,records,num_mm_f=0,num_mm_r=0):
    fp_df = pd.read_csv(fp_tab,sep='\t',header=None)
    rp_df = pd.read_csv(rp_tab,sep='\t',header=None)
    fp_df.columns = rp_df.columns = cols
    
    fp_dict = convert_df(fp_df)
    rp_dict = convert_df(rp_df)
    #print(fp_dict)
    all_products = []
    ## 
    fp_matches = filtering(fp_tab,fp_fasta,num_mismatch=num_mm_f)
    rp_matches = filtering(rp_tab,rp_fasta,num_mismatch=num_mm_r)
    matches = fp_matches.intersection(rp_matches)
#     print(matches,records)
    sub_r = [records[ti] for ti in matches]
#     print(sub_r)
    for r in sub_r:
        fp_s,fp_e = fp_dict[r.id]
        rp_s,rp_e = rp_dict[r.id]
        if fp_s <= rp_e:
            product = r[fp_s-1:rp_e]  # with primers
        else:
            product = r[rp_e-1:fp_s]  # with primers
        all_products.append((product,r.id))
    return all_products


def expand_dg_primer(seq):
   """return a list of all possible k-mers given a degenerate base"""
   return list(map("".join, product(*map(dgbase.get, seq))))

def write_degen_primers(target_seq,name_prefix,ofile,rc=False):
    if rc:
        target_seq = ''.join([IUPACData.ambiguous_dna_complement[_] for _ in target_seq[::-1]])
    a = expand_dg_primer(target_seq)
    with open(ofile,'w') as f:
        for i,seq in enumerate(a):
            f.write(f'>{name_prefix}\n{seq}\n')

#>nirSF
pf = """CCTAYTGGCCGCCRCART"""
#>nirSR
pr = """TCCMAGCCRCCRTCRTGCAG"""

reference_nirS = '/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/results/BACyqlbM/ref_nirS/correctAndRef_nirS_nucl.ffn'  
refernce_nirS_blastdb = './nirS.db'
fp_tab_target = 'Ruegeria_nirS_f_primer.blast'
rp_tab_target = 'Ruegeria_nirS_r_primer.blast'
os.system(f"makeblastdb -in {reference_nirS} -dbtype nucl -out {refernce_nirS_blastdb}")
####################
## write out degenated primer
pf_file = './nirS.F.fasta'
pr_file = './nirS.R.fasta'
write_degen_primers(pf,'nirSF',pf_file)
write_degen_primers(pr,'nirSR',pr_file)
####################
## blast against reference
cmds = f"""
blastn -task blastn-short  -query {pf_file} -db {refernce_nirS_blastdb} -outfmt '6 qseqid sseqid pident length evalue mismatch gapopen sstart send stitle' -num_threads 6  -out {fp_tab_target}
blastn -task blastn-short  -query {pr_file} -db {refernce_nirS_blastdb} -outfmt '6 qseqid sseqid pident length evalue mismatch gapopen sstart send stitle' -num_threads 6  -out {rp_tab_target}
"""
check_call(cmds,shell=1)
####################
## parse blast results and get amplified regions (product)
results = {}
records_target = {r.id:r for r in SeqIO.parse(reference_nirS,'fasta')}
for mm_pattern in [(0,0),(0,1),(1,0),(1,1)]:
    f,r = mm_pattern
    product_ = get_product(fp_tab_target,rp_tab_target,pf_file,pr_file,records_target,f,r)
    results[('mapped',f,r)] = product_
    print('target',mm_pattern,len(product_),len(records_target))
    
    
    