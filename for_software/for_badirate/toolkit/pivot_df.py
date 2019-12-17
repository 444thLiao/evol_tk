"""
Formatting results from kegg annotation(diamond output)
"""
from collections import defaultdict

import pandas as pd
from tqdm import tqdm

from dating_workflow.step_script import convert_genome_ID_rev

header = 'qseqid sseqid salltitles pident length evalue bitscore'.split(' ')

locusID2kegg_list = '/home-user/thliao/data/protein_db/kegg/latest/links/genes_ko.list'
locusID2kegg_dict = dict([_.split('\t') for _ in open(locusID2kegg_list).read().split('\n') if _])

infile = "./protein_annotations/kegg_merged.tab"
with open('./protein_annotations/kegg_merged_e50_top5.tab', 'w') as f1:
    remained_rows = []
    _count = 0
    used_locusID = set()
    for row in tqdm(open(infile)):
        rows = row.split('\t')
        row_dict = dict(zip(header, rows))
        if rows[0] in used_locusID:
            _count += 1
        else:
            _count = 0
            used_locusID.add(rows[0])
        if _count >= 50:
            continue
        if float(row_dict['evalue']) <= 1e-50:
            # remained_rows.append(row)
            f1.write(row)
            f1.flush()

used_dict = defaultdict(list)
infile = "./protein_annotations/kegg_merged.tab"
for row in tqdm(open(infile)):
    rows = row.split('\t')
    row_dict = dict(zip(header, rows))
    if len(used_dict[rows[0]]) > 10:
        continue
    if float(row_dict['evalue']) <= 1e-50:
        used_dict[rows[0]].append(locusID2kegg_dict.get(rows[1], None))

# all_ko = [_.split(':')[-1] for v in used_dict.values() for _ in v if _ is not None]
# all_ko = list(set(all_ko))

all_None_seqs = []
g2ko2tags = defaultdict(lambda: defaultdict(list))
for locus_tag, annotation in tqdm(used_dict.items()):
    gid = convert_genome_ID_rev(locus_tag)
    valid_annotations = list(set([_ for _ in annotation if _ is not None]))
    if len(valid_annotations) == 1:
        ko = valid_annotations[0].split(':')[-1]
        g2ko2tags[gid][ko].append(locus_tag)
    elif len(valid_annotations) > 1:
        for ko in set(valid_annotations):
            ko = ko.split(':')[0]
            g2ko2tags[gid][ko].append(locus_tag)
        # multi_match.append(locus_tag)
    elif len(valid_annotations) != 0:
        all_None_seqs.append(locus_tag)
    else:
        pass

g2ko2tags = {k: {_k: ','.join(_v) for _k, _v in v.items()}
             for k, v in g2ko2tags.items()}
total_df = pd.DataFrame.from_dict(g2ko2tags, orient='index')
# totaldf = total_df.fillna(0)
total_df.to_csv('./protein_annotations/kegg_diamond_info.crosstab', sep='\t', index=1)
