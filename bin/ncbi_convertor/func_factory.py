import io
import os
from os.path import exists, dirname

import click
from Bio import Entrez
from tqdm import tqdm

from bin.ncbi_convertor import edl, access_intermedia, parse_id



def aid2GI(id_list, redo=False):
    "prototype of function which convert ID to GI"
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
            _results, failed = edl.esearch(db='protein',
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
