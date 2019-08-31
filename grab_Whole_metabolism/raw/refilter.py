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

locus_list = [_.id for _ in collect_reads]
locus2gene_df = pd.read_csv("../N-relative_genes.tsv", sep='\t', index_col=0)
sub_df = pre_df.loc[pre_df.loc[:, 0].isin(locus_list), :]

locus2ko = dict()
locus2module = dict()
locus2completeOrthos = dict()

choose_highest_one = sub_df.sort_values([0, 10]).drop_duplicates(0)
for rid, row in tqdm(choose_highest_one.iterrows(),
                     total=choose_highest_one.shape[0]):
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
    if not isinstance(ko, str):
        ko = ';'.join(set(list(ko.values)))
    _sub_df = locus2gene_df.loc[locus2gene_df.loc[:, 'Orthology(single)'] == ko, :]
    name = [_
            for _ in _sub_df.Name.value_counts().index
            if _ != 'nan'][0]
    locus2name[locus] = name
name_counts = Counter([_ for _ in locus2name.values()])
############################################################
sample2info = pd.read_csv('sample2infos.tsv', sep='\t', header=0, index_col=1)

locus2info_df = pd.DataFrame(columns=["locus",
                                      "locus_prefix",
                                      'sample name',
                                      'source project',
                                      'ko(single)',
                                      'Gene name(N metabolism)',
                                      'ko(complete)',
                                      'module',
                                      ])
from collections import defaultdict

count_ = 0
locus2info_dict = defaultdict(dict)
for locus, g_name in tqdm(locus2name.items(),
                          total=len(locus2name)):
    locus_prefix = locus.split('_')[0]
    sname, source = sample2info.loc[locus_prefix, :].values
    if not isinstance(locus2ko[locus], str):
        for rid, v in enumerate(locus2ko[locus]):
            locus2info_df.loc[count_, :] = [locus,
                                            locus_prefix,
                                            sname,
                                            source,
                                            locus2ko[locus].values[rid],
                                            g_name,
                                            locus2completeOrthos[locus].values[rid],
                                            locus2module[locus].values[rid]]
            count_ += 1
    else:
        locus2info_df.loc[count_, :] = [locus,
                                        locus_prefix,
                                        sname,
                                        source,
                                        locus2ko[locus],
                                        g_name,
                                        locus2completeOrthos[locus],
                                        locus2module[locus]]
        count_ += 1
# manually curated
locus2info_df.loc[locus2info_df.loc[:, "Gene name(N metabolism)"] == 'NEUTE1DRAFT_87025', 'Gene name(N metabolism)'] = 'nit-6'
# for ko = K17877
locus2info_df.to_csv('./contain_N_relative_locus2info.tsv', sep='\t', index=0)

############################################################
from os.path import join, exists, dirname, basename
from subprocess import check_call
import os
from glob import glob
from tqdm import tqdm
import multiprocessing as mp


def run_cmd(cmd):
    check_call(cmd, shell=True)


cmd_template = "/home-user/thliao/bin/kraken2 --quick --db /home-backup/thliao/kraken2_db/k2db --threads 40 --report {outfile} --memory-mapping {infile} --output -"

base_dir = '/home-user/thliao/data/metagenomes/concat_all/'
for input_fna in tqdm(glob(join(base_dir, 'prokka_o', '*', '*.fna'))):
    g = basename(dirname(input_fna))
    os.makedirs(join(base_dir, 'k_output'), exist_ok=True)
    # input_fna = glob(join(base_dir, 'prokka_o', g, '*.fna'))[0]
    ofile = join(base_dir, 'k_output', g + '.kout')
    if not exists(ofile):
        run_cmd(cmd_template.format(infile=input_fna,
                                    outfile=ofile))
kraken2_header = ["percentage_frag",
                  "num frag",
                  "num assigned frag",
                  "rank code",
                  "NCBI taxid",
                  "scientific name"]
levels = ["D", "P", "C", "O", "F", "S"]


def parse_kraken2(infile):
    # todo: check some abnormal situations.
    df = pd.read_csv(infile, sep='\t', header=None)
    df.columns = kraken2_header
    df.loc[:, "scientific name"] = [_.strip()
                                    for _ in df.loc[:, "scientific name"]]
    sorted_df = df.sort_values("percentage_frag", ascending=False)
    # sorted_df = sorted_df.loc[df.loc[:, "rank code"] == "S", :]
    for l in levels[::-1]:
        _df = sorted_df.loc[sorted_df.loc[:, 'rank code'] == l, :]
        if not _df.shape[0]:
            continue
        if _df.iloc[0, 0] <= 20:
            continue
        return _df


from ete3 import NCBITaxa

ncbi = NCBITaxa()

sample2infos = pd.read_csv('/home-user/thliao/data/metagenomes/concat_all/sample2infos.tsv', sep='\t', index_col=0)
sample2infos = sample2infos.reindex(columns=['locus_prefix',
                                             'source',
                                             'superkingdom',
                                             'phylum',
                                             'class',
                                             'order',
                                             'family',
                                             "genus",
                                             'species',
                                             'superkingdom(from metadata)',
                                             'phylum(from metadata)',
                                             'class(from metadata)',
                                             'order(from metadata)',
                                             'family(from metadata)',
                                             'genus(from metadata)',
                                             'species(from metadata)', ])

for g in tqdm(sample2infos.index):
    ofile = join(base_dir, 'k_output', g + '.kout')
    df = parse_kraken2(ofile)
    if df is None:
        continue
    tid = df.iloc[0, 4]
    lineage = ncbi.get_lineage(tid)
    rank = ncbi.get_rank(lineage)
    rank = {v: k for k, v in rank.items()}
    names = ncbi.get_taxid_translator(lineage)
    for c in ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
        if c in rank:
            sample2infos.loc[g, c] = names[rank[c]]

s2df = {}
for source in sample2infos.source.unique():
    target_file = f'/home-user/thliao/data/metagenomes/{source}/metadata.csv'
    if exists(target_file):
        metadata = pd.read_csv(target_file, sep='\t')
    else:
        metadata = None
    s2df[source] = metadata
no_metadata = """19_Stewart
17_lee
18_Delmont
"""

tid_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
for source, df in s2df.items():
    if df is not None:
        for _, row in tqdm(df.iterrows(), total=df.shape[0]):
            sid = row['assembly_accession']
            tid = row['taxid']
            lineage = ncbi.get_lineage(tid)
            rank = ncbi.get_rank(lineage)
            rank = {v: k for k, v in rank.items()}
            names = ncbi.get_taxid_translator(lineage)
            for tlevel in tid_levels:
                if tlevel in rank:
                    sample2infos.loc[sample2infos.index.str.startswith(sid),
                                     tlevel + '(from metadata)'] = names[rank[tlevel]]

##  manually assigned......tired
t = pd.read_excel('19_Stewart/metadata.xlsx')
tids = [row['original_bin']
        for rid, row in tqdm(t.iterrows())]
t = t.set_index('original_bin')
t = t.iloc[:,[11,12,13,14,15,16,17]]
t = pd.DataFrame(t.values,
                 columns= ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'],
                 index=t.index)
print(sample2infos.index[sample2infos.source == '19_Stewart'].difference(set(tids)))
tids = sample2infos.index[sample2infos.source == '19_Stewart'].intersection(set(tids))

sample2infos.loc[tids,
                 [tlevel + '(from metadata)'
                  for tlevel in tid_levels]] = t.loc[tids, tid_levels].values
##
t = pd.read_excel('17_lee/metadata.xlsx', sheet_name=1,header=1)
t = t.iloc[:, [0, 17, 18, 19, 20, 21, 22, 23]]
t.columns = ['bin_id', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
t = t.set_index('bin_id')
tids = [_+'-contigs' for _ in list(t.index)]
t.index = tids
print(sample2infos.index[sample2infos.source == '17_lee'].difference(set(tids)))
tids = sample2infos.index[sample2infos.source == '17_lee'].intersection(set(tids))

sample2infos.loc[tids,
                 [tlevel + '(from metadata)'
                  for tlevel in tid_levels]] = t.loc[tids, tid_levels].values
##
t = pd.read_excel('18_Delmont/MAG_metadata(new).xlsx', sheet_name=0)
t = t.iloc[:, [0, 15, 16, 17, 18, 19, 20, 21]]
t.columns = ['bin_id', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
t = t.set_index('bin_id')
tids= list(t.index)
print(sample2infos.index[sample2infos.source == '18_Delmont'].difference(set(tids)))
tids = sample2infos.index[sample2infos.source == '18_Delmont'].intersection(set(tids))
sample2infos.loc[tids,
                 [tlevel + '(from metadata)'
                  for tlevel in tid_levels]] = t.loc[tids, tid_levels].values

############################################################
# summary
#

# GET a module based df
# from collections import Counter
#
# t = pd.read_csv("N-relative_genes.tsv", sep='\t')
# for m in t.loc[:, 'module Name'].unique():
#     sub_t = t.loc[t.loc[:, 'module Name'] == m, :]
#     print(m, Counter(sub_t.Name.fillna(0)))


# def classified(sub_df):
#     # for a tax level
#     idx2module = {0: 'Dissimilatory nitrate reduction, nitrate => ammonia',
#                   1: 'Denitrification, nitrate => nitrogen',
#                   2: 'Complete nitrification, comammox, ammonia => nitrite => nitrate',
#                   3: 'Assimilatory nitrate reduction, nitrate => ammonia',
#                   4: 'Nitrogen fixation, nitrogen => ammonia',
#                   5: 'Nitrification, ammonia => nitrite'}
#
#     genes = sub_df.loc[:, "Gene name(N metabolism)"]
#     idx2genes = {4:{'nifH','nifD','nifK',
#                     'anfH','anfD','anfK','anfG',
#                     'vnfH','vnfD','vnfK','vnfG'},
#                  0:{'narH','narG',
#                     'nasA','nasB',
#                     'nirA',
#                     'nrfA'}}

classified = [["Nitrogen Fixation", "nifD", "nifK", "nifH", "anfG"],
              ["Assimilatory Nitrate Reduction", "narB", "NR", "nasA", "nasB"],
              ["Assimilatory Nitrite Reduction", "nit-6", "nirA"],
              ["Dissimilatory Nitrate Reduction", "narG", "narH", "narI", "napA", "napB"],
              ["Dissimilatory Nitrite Reduction", "nirB", "nirD", "nrfA", "nrfH"],
              ["Denitrification", "nirK", "nirS", "norB", "norC", "norH", "nosZ"],
              ["Partial Denitrification (NO2- to NO)", "nirK", "nirS"],
              ["Partial Denitrification (NO to N2O)", "norB", "norC"],
              ["Partial Denitrification (N2O to N2)", "nosZ"],
              ["Nitrification", "amoC", "amoA", "amoB", "hao"],
              ["Ammonium to Hydroxylamine", "amoC", "amoA", "amoB"]]

# bar plot for phylum
import plotly
from plotly import graph_objs as go

level = 'phylum'

fig = go.Figure()
sub_df = locus2info_df.copy()
sub_df.loc[:, level].replace('', 'unclassified', inplace=True)
for m in sub_df.loc[:, 'module'].unique():
    _df = sub_df.loc[sub_df.module == m, :]
    count_data = _df.loc[:, level].value_counts()
    total_count = sub_df.loc[:, level].value_counts()
    freq_data = count_data / total_count * 100
    fig.add_trace(go.Bar(x=count_data.index,
                         y=count_data.values,
                         name=m))
