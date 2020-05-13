"""
It mainly compares the differences between annotations against different database.

"""
from tqdm import tqdm
from os.path import expanduser
import pandas as pd
kola = open(expanduser('~/Rhizobiales/kegg_GhostKOALA/annotated.txt'))
kola_dict = {}
for row in tqdm(kola):
    rows = row.strip('\n').split('\t')
    if len(rows)==2:
        kola_dict[rows[0]] = rows[1]
    elif len(rows) > 2:
        print(rows)

from glob import glob

header = 'qseqid sseqid pident length evalue bitscore'.split(' ')

locusID2kegg_list = '/home-user/thliao/data/protein_db/kegg/latest/links/genes_ko.list'
locusID2kegg_dict = dict([_.split('\t') for _ in open(locusID2kegg_list).read().split('\n') if _])

total_count = 0
fall = open('/home-user/jjtao/Rhizobiales/diamond-test/kegg_diamond_merged_e50_top50.tab', 'w')
for f in tqdm(glob("/home-user/jjtao/Rhizobiales/diamond-test/kegg_diamond/*.blast8")):
    _count = 0
    used_locusID = set()
    for row in tqdm(open(f)):
        rows = row.split('\t')
        row_dict = dict(zip(header, rows))
        if rows[0] in used_locusID:
            _count += 1
        else:
            _count = 0
            used_locusID.add(rows[0])
        if _count >= 50:
            continue
        if float(rows[-2]) <= 1e-50:
            total_count +=1
            fall.write(row)
            fall.flush()

from collections import defaultdict
used_dict = defaultdict(list)
infile = "/home-user/jjtao/Rhizobiales/diamond-test/kegg_diamond_merged_e50_top50.tab"
for row in tqdm(open(infile)):
    rows = row.split('\t')
    row_dict = dict(zip(header, rows))
    used_dict[rows[0]].append((locusID2kegg_dict.get(rows[1], None),float(rows[-2])))

tag2ko = {}
for locus_tag, annotation in tqdm(used_dict.items()):
    valid_annotations = list(set([_ for _ in annotation if _[0] is not None]))
    if len(valid_annotations) == 1:
        ko = valid_annotations[0]
        tag2ko[locus_tag] = ko[0].split(':')[-1]
    elif len(valid_annotations) > 1:
        smallest = 100
        smallest_ko = ''
        for ko,e in set(valid_annotations):
            if e <= smallest:
                ko = ko.split(':')[0]
                smallest_ko = ko
                smallest = e
            tag2ko[locus_tag] = smallest_ko




dup_ids = open('/home-user/jjtao/Rhizobiales/kegg_hmmsearch/dup_ids.list').read().split('\n')
dup_ids = set(dup_ids)
hmm_dict = {}
_df3 = pd.read_csv('/home-user/jjtao/Rhizobiales/kegg_hmmsearch/merged_result/merged_hmm_info.tab',sep='\t',index_col=0)
for ko,row in tqdm(_df3.iterrows()):
    for _ in row:
        if isinstance(_,str):
            for v in _.split(','):
                if v in dup_ids:
                    v = v.rpartition('|')[0]
                hmm_dict[v] = ko

df1 = pd.DataFrame().from_dict(kola_dict,orient='index')
df2 = pd.DataFrame().from_dict(tag2ko,orient='index')
df3 = pd.DataFrame().from_dict(hmm_dict,orient='index')
df1.columns = ['kola_KO']
df2.columns = ['blast_KO']
df3.columns = ['hmm_KO']
merged_df = pd.concat([df1,df2,df3],axis=1)

#############
def compare_two(a,b,total=2563959):
    col1 = merged_df.loc[:, a].dropna()
    col2 = merged_df.loc[:, b].dropna()
    # total_annotated
    num_col1 = round(len(col1)/total*100,2)
    num_col2 = round(len(col2)/total*100,2)
    # num_shared
    num_shared = round(len(col1.index.intersection(col2.index))/total*100,2)
    # identital
    num_ident = round((col2.reindex(col1.index) == col1).sum()/total*100,2)

    # num col1 but not col2
    at_col1_not_at_col2 = round((num_col1 - num_shared),2)
    at_col2_not_at_col1 = round((num_col2 - num_shared),2)

    print(f"""
{a} annotated {num_col1}% genes
{b} annotated {num_col2}% genes
both {a}{b} annotated {num_shared}% genes
Among above both annotated genes
{num_ident}% (versus all genes instead of both annotated) genes are identitcal
{at_col1_not_at_col2}% genes annotated by {a} but not annotated by {b}
{at_col2_not_at_col1}% genes annotated by {b} but not annotated by {a}
""")

import itertools
for a,b in list(itertools.combinations(merged_df.columns,2)):
    print('')
    print('')
    print('')
    compare_two(a,b)