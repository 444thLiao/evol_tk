from collections import defaultdict
from ete3 import NCBITaxa
from Bio import Entrez
from tqdm import tqdm

ncbi = NCBITaxa()
Entrez.email = 'l0404th@gmail.com'
# taxids = 1224 # Proteobacteria

all_levels_collect = defaultdict(list)
taxons = ['phylum', 'class',
          'order', 'family', 'genus']

name = ''
names = [_.strip() for v in name.split('|') for _ in v.split(',')]
def descendent_tax():
    lineage = ncbi.get_des

def name2aid(tname):
    try:
        result = Entrez.read(Entrez.esearch(db='assembly', term=tname, retmax=100, sort='Significance'))
        get_id = result['IdList'][0]
        assembly_id = Entrez.read(Entrez.esummary(db='assembly', id=get_id))
        aid = assembly_id['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Genbank']
    except:
        print(tname)
        return
    #aid = aid.replace('GCF', 'GCA')
    return aid

for tid in taxids:
    lineage = ncbi.get_lineage(tid)
    ranks = ncbi.get_rank(lineage)
    names = ncbi.get_taxid_translator(lineage)
    used_ranks = [_t for _t, r in ranks.items() if r in taxons]
    for r in used_ranks:
        rank = ranks[r]
        all_levels_collect[rank].append(names[r])
all_levels_collect = {k: set(v) for k, v in all_levels_collect.items()}

collect_tid = []
for p, v in tqdm(all_levels_collect.items()):
    #
    for tname in v:
        result = Entrez.read(Entrez.esearch(db='assembly', term=tname, retmax=100, sort='Significance'))
        get_id = result['IdList'][0]
        assembly_id = Entrez.read(Entrez.esummary(db='assembly', id=get_id))
        aid = assembly_id['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']
        aid = aid.replace('GCF', 'GCA')
        collect_tid.append(aid)
