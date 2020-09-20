"""
RAW scripts , tries

"""

from ForOrthofinder.bin.resort_OG_with_gaps import *

infile = "./sorted_OG_renamed.csv"
backbone_column_idx = 0

OG_df = pd.read_csv(infile,sep='\t',index_col=0)
genome2order_tuple = pickle.load(open(join(tmp_dir, 'genome2order_tuple'), 'rb'))
genome2gene_info = pickle.load(open(join(tmp_dir, 'genome2gene_info'), 'rb'))

# raw part (rename)
raw_df = pd.read_csv("./cat_sub_Orthogroups.csv",sep='\t',index_col=0)
locus_id2name = {}
for col in tqdm(raw_df.columns):
    col_v = raw_df[col]
    col_v = [_ for _ in col_v if not pd.isna(_)]
    for _ in col_v:
        vals = _.split(',')
        for val in vals:
            locus_id = val.split('|')[-1].split(' ')[0].strip()
            locus_id2name[locus_id] = val
def f(x):
    if pd.isna(x):
        return x
    else:
        locus_id = x.split('|')[-1].split(' ')[0].strip()
        if locus_id2name[locus_id] != x:
            return locus_id2name[locus_id]
        else:
            return x
new_OG_df = OG_df.applymap(f)
# raw part
def get_contig(x,gene_info):
    # unique to luo lab faa format
    locus = x.split('|')[-1].split(' ')[0].strip()
    info = gene_info[locus]
    return info['contig_name']

backbone_col = OG_df.columns[backbone_column_idx]
resorted_rows = OG_df.loc[OG_df.loc[:,backbone_col].isna(),:]
OG_df.loc[:,'reliability'] = ''
OG_df.loc[~OG_df.loc[:,backbone_col].isna(),'reliability'] = True

reliability = []
for idx, row in tqdm(resorted_rows.iterrows()):
    used_genome = row.index[~row.isna()][0]
    v = row[used_genome]

    num_idx = OG_df.index.get_loc(idx)
    up_right = ''
    down_right = ''
    for _ in OG_df.index[:num_idx][::-1]:
        if _ not in resorted_rows.index and not pd.isna(OG_df.loc[_, used_genome]):
            up_right = _
            break
    for _ in OG_df.index[num_idx + 1:]:
        if _ not in resorted_rows.index and not pd.isna(OG_df.loc[_, used_genome]):
            down_right = _
            break
    # up_right = [_ for _ in OG_df.index[:num_idx] if _ not in resorted_rows.index and not pd.isna(OG_df.loc[_,used_genome])][-1]
    # down_right = [_ for _ in OG_df.index[num_idx+1:] if _ not in resorted_rows.index and not pd.isna(OG_df.loc[_,used_genome])][0]
    v_up = OG_df.loc[up_right, used_genome]
    v_down = OG_df.loc[down_right, used_genome]
    gene_info = genome2gene_info[used_genome]
    v_c = get_contig(v,gene_info)
    v_up_c = get_contig(v_up,gene_info)
    v_down_c = get_contig(v_down, gene_info)
    if len(set([v_c,v_up_c,v_down_c])) <=2:
        OG_df.loc[idx,'reliability'] = True
        reliability.append(True)
    else:
        OG_df.loc[idx, 'reliability'] = False
        reliability.append(False)
