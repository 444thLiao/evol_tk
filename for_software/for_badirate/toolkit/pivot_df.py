"""
Formatting results from kegg annotation(diamond output)
"""
import pandas as pd
from collections import defaultdict
from tqdm import tqdm
from dating_workflow import convert_genome_ID,convert_genome_ID_rev   

header = 'qseqid sseqid salltitles pident length evalue bitscore'.split(' ')

locusID2kegg_list = '/home-user/thliao/data/protein_db/kegg/latest/links/genes_ko.list'
locusID2kegg_dict = dict([_.split('\t') for _ in open(locusID2kegg_list).read().split('\n') if _ ] )

infile = "./protein_annotations/kegg_merged.tab"
with open('./protein_annotations/kegg_merged_filtered.tab','w') as f1:
    remained_rows = []
    for row in tqdm(open(infile)):
        rows = row.split('\t')
        row_dict = dict(zip(header,rows))
        if float(row_dict['evalue']) <= 1e-50:
            # remained_rows.append(row)
            f1.write(row)
            f1.flush()

used_dict = defaultdict(list)
infile = "./protein_annotations/kegg_merged.tab"
for row in tqdm(open(infile)):
    rows = row.split('\t')
    row_dict = dict(zip(header,rows))
    if len(used_dict[rows[0]])>5:
        continue
    if float(row_dict['evalue']) <= 1e-50:
        used_dict[rows[0]].append(locusID2kegg_dict.get(rows[1]))

# all_ko = [_.split(':')[-1] for v in used_dict.values() for _ in v if _ is not None]
# all_ko = list(set(all_ko))


multi_match = []
g2ko2num = defaultdict(lambda : defaultdict(int))
for locus_tag,annotation in tqdm(used_dict.items()):
    gid = convert_genome_ID_rev(locus_tag)
    valid_annotations = list(set([_ for _ in annotation if _ is not None]))
    if len(valid_annotations) == 1:
        ko = valid_annotations[0].split(':')[-1]
        g2ko2num[gid][ko] += 1
    elif len(valid_annotations) > 1:
        multi_match.append(locus_tag)
    else:
        pass
    
total_df = pd.DataFrame.from_dict(g2ko2num,orient='index')
totaldf = total_df.fillna(0)
totaldf.to_csv('./protein_annotations/kegg_diamond.crosstab',sep='\t',index=1)



