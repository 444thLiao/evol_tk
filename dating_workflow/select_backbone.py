from collections import defaultdict

from Bio import Entrez
from ete3 import NCBITaxa
from tqdm import tqdm

ncbi = NCBITaxa()
Entrez.email = 'l0404th@gmail.com'
# taxids = 1224 # Proteobacteria


taxons = ['phylum', 'class',
          'order', 'family', 'genus', 'species']

name = ''
names = [_.strip() for v in name.split('|') for _ in v.split(',')]


def get_descendent_tax(tid, next_level):
    children_ids = ncbi.get_descendant_taxa(tid)
    collect_tid = []
    for cid in children_ids:
        lineage = ncbi.get_lineage(cid)
        levels_dict = ncbi.get_rank(lineage)
        levels_dict = {v: k for k, v in levels_dict.items()}
        _cache = levels_dict.get(next_level, '')
        if _cache:
            collect_tid.append(_cache)
    return set(collect_tid)


def name2aid(tname):
    try:
        result = Entrez.read(Entrez.esearch(db='assembly', term=tname, retmax=100, sort='Significance'))
        get_id = result['IdList'][0]
        assembly_id = Entrez.read(Entrez.esummary(db='assembly', id=get_id))
        aid = assembly_id['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Genbank']
    except:
        # pass
        # print(tname)
        return
    # aid = aid.replace('GCF', 'GCA')
    return aid


taxids = []
for row in open('./concat_metadata.csv'):
    if row:
        taxids.append(row.split('\t')[6])
taxids = taxids[1:]
taxids = set(taxids)
all_levels_collect = defaultdict(list)
for tid in taxids:
    lineage = ncbi.get_lineage(tid)
    ranks = ncbi.get_rank(lineage)
    names = ncbi.get_taxid_translator(lineage)
    used_ranks = [_t for _t, r in ranks.items() if r in taxons[:-1]]
    for r in used_ranks:
        rank = ranks[r]
        all_levels_collect[rank].append((names[r], r))
all_levels_collect = {k: set(v) for k, v in all_levels_collect.items()}

prepared_tids = []
for level, tids in tqdm(all_levels_collect.items()):
    next_level = taxons[taxons.index(level) + 1]
    tmp_collect = set()
    for tid in tqdm(tids):
        next_tids = get_descendent_tax(int(tid[1]), next_level)
        tmp_collect.update(next_tids)
    prepared_tids += list(tmp_collect)

prepared_tids = set(prepared_tids)
# prepared_tids to representative genomes
tid2aid = {}
tnames = ncbi.get_taxid_translator(prepared_tids)
for tid, tname in tqdm(tnames.items()):
    aid = name2aid(tname)
    if aid is not None:
        tid2aid[tid] = aid
