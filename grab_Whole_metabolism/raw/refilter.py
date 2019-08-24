import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
from Bio import SeqIO
from collections import Counter

whole_kegg = 'N_relative_2_whole_kegg_blastp.out'
N_relative = 'N_relative_blastp.out'

pre_df = pd.read_csv(N_relative, sep='\t', header=None)
aft_df = pd.read_csv(whole_kegg, sep='\t', header=None)


def cal_ratio(locus):
    setA = set(pre_df.loc[pre_df.loc[:, 0] == locus, 1])
    setB = set(aft_df.loc[aft_df.loc[:, 0] == locus, 1])
    intersec = setA.intersection(setB)
    union_set = setA.union(setB)
    if len(intersec) / len(union_set) >= 0.5:
        # real_N_metabolism_genes.append(locus)
        return locus


real_N_metabolism_genes = []
locus_set = list(set(pre_df.loc[:, 0]))
with mp.Pool(processes=50) as tp:
    for locus in tqdm(tp.imap(cal_ratio, locus_set),
                      total=len(locus_set)):
        if locus is not None:
            real_N_metabolism_genes.append(locus)

records = SeqIO.parse('N_relative_blastp.faa', format='fasta')
collect_reads = [_ for _ in records if _.id in real_N_metabolism_genes]

with open('real_N_metabolism.faa', 'w') as f1:
    SeqIO.write(collect_reads, f1, format='fasta-2line')

locus2gene_df = pd.read_csv("../N-relative_genes.csv", sep='\t', index_col=0)
sub_df = pre_df.loc[pre_df.loc[:, 0].isin(locus_list), :]

locus2ko = dict()
locus2module = dict()
locus2completeOrthos = dict()

choose_highest_one = sub_df.sort_values([0, 10]).drop_duplicates(0)
for rid, row in tqdm(choose_highest_one.iterrows()):
    locus_tag = row[1]
    seq_name = row[0]
    name = locus2gene_df.loc[locus_tag, "Name"]
    ortho = locus2gene_df.loc[locus_tag, "Orthology(single)"]
    module = locus2gene_df.loc[locus_tag, "module Name"]
    completeOrthos = locus2gene_df.loc[locus_tag, "Orthology(total)"]

    locus2ko[seq_name] = ortho
    locus2module[seq_name] = module
    locus2completeOrthos[seq_name] = completeOrthos
locus2name = {}
# module_counts = Counter([tuple(sorted(_)) for _ in locus2module.values()])
for locus, ko in tqdm(locus2ko.items()):
    _sub_df = locus2gene_df.loc[locus2gene_df.loc[:, 'Orthology(single)'] == ko, :]
    name = [_
            for _ in _sub_df.Name.value_counts().index
            if _ != 'nan'][0]
    locus2name[locus] = name
name_counts = Counter([_ for _ in locus2name.values()])
############################################################
sample2info = pd.read_csv('sample2infos.tsv',sep='\t',header=0,index_col=1)

from collections import defaultdict
locus2info_dict = defaultdict(dict)
for locus in locus2name.keys():
    locus_prefix = locus.split('_')[0]
    sname,source = sample2info.loc[locus_prefix,:].values

    locus2info_dict[locus]['locus_prefix'] = locus_prefix
    locus2info_dict[locus]['sample name'] = sname
    locus2info_dict[locus]['source project'] = source
    locus2info_dict[locus]['ko(single)'] = locus2ko[locus]
    locus2info_dict[locus]['Gene name(N metabolism)'] = locus2name[locus]
    locus2info_dict[locus]['ko(complete)'] = locus2completeOrthos[locus]
    locus2info_dict[locus]['module'] = locus2module[locus]
locus2info_df = pd.DataFrame.from_dict(locus2info_dict,orient='index')
locus2info_df.to_csv('./contain_N_relative_locus2info.tsv',sep='\t',index=1,index_label="locus_tag")