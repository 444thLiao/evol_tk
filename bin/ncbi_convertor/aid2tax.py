"""
Assembly id to taxnomy id

"""
import io
import os
from os.path import exists, dirname

import click
from Bio import Entrez
from tqdm import tqdm

from bin.ncbi_convertor import edl, access_intermedia, parse_id
from ete3 import ncbi


def aid2GI(id_list, redo=False):
    suffix = 'aid2GI'
    _cache = access_intermedia(id_list, suffix=suffix, redo=redo)
    if _cache is not None:
        id2gi = _cache
    else:
        tqdm.write('from protein Accession ID to GI')
        results, failed = edl.esearch(db='assembly',
                                      ids=id_list,
                                      result_func=lambda x: Entrez.read(
                                          io.StringIO(x))['IdList'],
                                      batch_size=1
                                      )
        _results_dict = {}
        _count = 0
        while failed:
            failed_id_list = failed
            _results, failed = edl.esearch(db='assembly',
                                           ids=failed_id_list,
                                           result_func=lambda x: Entrez.read(
                                               io.StringIO(x))['IdList'],
                                           batch_size=1
                                           )
            _results_dict.update(dict(_results))
            _count += 1
            if _count >= 5:
                break
        tqdm.write('still %s failed IDs, be careful.....' % len(failed))
        # for edl.esearch, it will auto **zip** searched term and its result.
        id2gi = dict(results)
        id2gi.update(_results_dict)
        id2gi = {pid: id2gi.get(pid, '') for pid in id_list}
        # stodge the result into intermedia file for second access.
        access_intermedia(id2gi, suffix=suffix)
    return id2gi



id2gi = dict(results)
all_GI = list(set(id2gi.values()))

results, failed = edl.esummary(db='assembly',
                                ids=all_GI,)
new_results = []
for r in results:
    d = Entrez.read(io.StringIO(r),validate=False)
    new_results.extend(d['DocumentSummarySet']['DocumentSummary'])
    
aid2sci_name = {}
for r in new_results:
    aid = r['AssemblyAccession'].replace('GCF','GCA')
    name = r['Organism']
    aid2sci_name[aid] = name
    
    
aid2tid = {}
for r in new_results:
    aid = r['AssemblyAccession'].replace('GCF','GCA')
    tid = r['Taxid']
    aid2tid[aid] = tid
    
    lineage = ncbi.get_lineage(int(tid))
    rank = ncbi.get_rank(lineage)
    rank = {v: k for k, v in rank.items()}
    names = ncbi.get_taxid_translator(lineage)
    _dict = {}
    for c in taxons:
        if c in rank:
            _dict[c] = names[rank[c]]
    genome2tax.update( {aid.split('.')[0]: _dict} )
for c in taxons:
    if c in rank:
        pid2info_dict[each_aid][c] = names[rank[c]]









