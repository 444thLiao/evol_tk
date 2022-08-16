
import pandas as pd
from Bio import SeqIO

cols = ['qaccver',
 'saccver',
 'pident',
 'length',
 'mismatch',
 'gapopen',
 'qstart',
 'qend',
 'sstart',
 'send',
 'evalue',
 'bitscore']

fp_fasta = 'forward_degen.fasta'
rp_fasta = 'reverse_degen.fasta'

fp_tab = 'out_forward.txt'
rp_tab = 'out_reverse.txt'

ref_fasta = '/home-user/thliao/db/rpoB/matched_rpoB_nucleotide.fasta'

# simple version
# def filtering(tab_path,iden=100):
#     df = pd.read_csv(tab_path,sep='\t',header=None)
#     sub_df = df[1][df[2]>=iden]
#     return set(sub_df)

# complex version
def filtering(tab_path,fa_file,num_mismatch=1):
    query_len = len(list(SeqIO.parse(fa_file,'fasta'))[0].seq)
    df = pd.read_csv(tab_path,sep='\t',header=None)
    df.columns = cols
    df.loc[:,'nm'] = (query_len - df['length']*df['pident']/100)
    sub_df = df['saccver'][df['nm']<=num_mismatch]
    return set(sub_df)

def convert_df(df):
    _tmp = df[['saccver','sstart','send']]
    _tmp = _tmp.set_index('saccver')
    _tmp = _tmp.loc[~_tmp.index.duplicated(),:]
    d = _tmp.to_dict('index')
    d = {k: (min([v['sstart'],v['send']]),
             max([v['sstart'],v['send']]) )  
         for k,v in d.items()}
    return d

def get_product(fp_tab,rp_tab,target_ids,records):
    fp_df = pd.read_csv(fp_tab,sep='\t',header=None)
    rp_df = pd.read_csv(rp_tab,sep='\t',header=None)
    fp_df.columns = rp_df.columns = cols
    
    fp_dict = convert_df(fp_df)
    rp_dict = convert_df(rp_df)
    
    all_products = []
    ## 
    sub_r = [records[ti] for ti in target_ids]
    for r in sub_r:
        fp_s,fp_e = fp_dict[r.id]
        rp_s,rp_e = rp_dict[r.id]
        if fp_s <= rp_e:
            product = r[fp_s-1:rp_e]  # with primers
        else:
            product = r[rp_e-1:fp_s]  # with primers
        all_products.append(product)
    return all_products
    
fp_matches = filtering(fp_tab,fp_fasta,num_mismatch=1)
rp_matches = filtering(rp_tab,rp_fasta,num_mismatch=1)

matches = fp_matches.intersection(rp_matches)
records = {r.id:r for r in SeqIO.parse(ref_fasta,'fasta')}

product = get_product(fp_tab,rp_tab,matches,records)

# show sequence description
for _ in product:
    print(' '.join(_.description.split(' ')[2:]))

# with open('./retrieved_products.fas','w') as f1:
#     SeqIO.write(product,f1,'fasta-2line')




