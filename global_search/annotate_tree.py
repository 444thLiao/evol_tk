from ete3 import Tree,NCBITaxa
# from global_search.thirty_party.EntrezDownloader import EntrezDownloader
from Bio import Entrez
import io
from tqdm import tqdm
from os.path import *
from subprocess import check_call
ncbi = NCBITaxa()

# edl = EntrezDownloader(
#     # An email address. You might get blocked by the NCBI without specifying one.
#     email='l0404th@gmail.com',
#     # An API key. You can obtain one by creating an NCBI account. Speeds things up.
#     api_key='ccf9847611deebe1446b9814a356f14cde08',
#     num_threads=30,                       # The number of parallel requests to make
#     # The number of IDs to fetch per request
#     batch_size=500,
#     pbar=True                             # Enables a progress bar, requires tqdm package
# )
    
basedir = '/home-user/thliao/project/nitrogen_cycle/fetch_genes/query_result'

tree_file = join(basedir,'nr_retrieve_hao/filtered_by_kegg.faa_aln.dir/iqtree.no_trim.treefile/K10535.sorted.newick')

def reformat(s):
    a = s.split('_')[-1]
    if not '_' in s:
        return s
    try:
        float(a)
        return s
    except:
        return s.rpartition('_')[0]
t = Tree(tree_file,format=1)
all_ids = t.get_leaf_names()
all_ids = [reformat(_) for _ in all_ids]
basedir = dirname(tree_file)
with open(join(basedir,'used_ids.list'),'w') as f1:
    f1.write('\n'.join(all_ids))

infile = join(basedir,'used_ids.list')
check_call(f'python3 /home-user/thliao/script/evolution_relative/global_search/get_info.py -i {infile} -o {dirname(infile)}',shell=1)

    
# results, failed = edl.esearch(db='protein',
#                                 ids=all_ids,
#                                 result_func=lambda x: Entrez.read(io.StringIO(x))['IdList'])
# all_GI = list(set(results[::]))
# tqdm.write('get pid summary from each one')
# results, failed = edl.esummary(db='protein',
#                                 ids=all_GI,
#                                 result_func=lambda x: Entrez.read(
#                                     io.StringIO(x)))
# id2complete_tax = {}
# id2tax = {}
# id2org = {}
# for r in tqdm(results):
#     aid = r['AccessionVersion']
#     tid = r['TaxId'].real
#     lineage = ncbi.get_lineage(tid)
#     rank = ncbi.get_rank(lineage)
#     rank = {v: k for k, v in rank.items()}
#     names = ncbi.get_taxid_translator(lineage)
#     rank2names = {r:names.get(tax,'') for r,tax in rank.items() if r !='no rank'}
#     id2complete_tax[aid] = rank2names
    
#     if names.get(rank.get('phylum', ''), 'ENV') == 'Proteobacteria':
#         id2tax[aid] = names.get(rank.get('class', ''), 'ENV')
#     else:
#         id2tax[aid] = names.get(rank.get('phylum', ''), 'ENV')
#     id2org[aid] = names[tid]
# id2tax = {k: v for k, v in id2tax.items() if k in final_ids}