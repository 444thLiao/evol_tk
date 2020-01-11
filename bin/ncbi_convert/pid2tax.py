# -*- coding: utf-8 -*-
"""
This script is mainly for implementing taxonomy information with its protein accession id.
It will generate a tab separated table.
"""
import io
import os
from collections import defaultdict
from os.path import exists, dirname

import click
from Bio import Entrez
from ete3 import NCBITaxa
from tqdm import tqdm

from bin.ncbi_convert import edl, access_intermedia, parse_id, taxons
from bin.ncbi_convert.pid2GI import pid2GI

ncbi = NCBITaxa()


def GI2tax(id2gi, redo=False):
    suffix = 'pid2tax'
    id_list = list(id2gi)
    _cache = access_intermedia(id_list, suffix=suffix, redo=redo)
    _results = []
    if _cache is not None:
        pid2info_dict = _cache
        return pid2info_dict
    all_GI = list(set(id2gi.values()))
    tqdm.write('get pid summary from each one')
    results, failed = edl.esummary(db='protein',
                                   ids=all_GI,
                                   result_func=lambda x: Entrez.read(
                                       io.StringIO(x)))
    if failed:
        tqdm.write("failed retrieve summary of %s protein ID" % len(failed))
        tqdm.write("retrieve each failed GI one by one")
        _results, _failed = edl.esummary(db='protein',
                                         ids=failed,
                                         batch_size=1,
                                         result_func=lambda x: Entrez.read(
                                             io.StringIO(x)))
        tqdm.write("failed ID is " + '\n'.join(map(str, _failed)))
    tqdm.write('from summary to GI and taxonomy')
    pid2info_dict = defaultdict(dict)
    for result in tqdm(results + _results):
        aid = result['AccessionVersion']
        _aid = [_ for _ in id2gi if _ in aid]
        # sometime, pid maybe a part of some intace one.
        if not _aid:
            continue
        else:
            right_aids = _aid
        for each_aid in right_aids:
            pid2info_dict[each_aid]['GI'] = result['Gi'].real
            pid2info_dict[each_aid]['taxid'] = result['TaxId'].real
            pid2info_dict[each_aid]['accession'] = aid
            try:
                lineage = ncbi.get_lineage(result['TaxId'].real)
                rank = ncbi.get_rank(lineage)
                rank = {v: k for k, v in rank.items()}
                names = ncbi.get_taxid_translator(lineage)
                for c in taxons:
                    if c in rank:
                        pid2info_dict[each_aid][c] = names[rank[c]]
            except:
                tqdm.write("failed to parse taxonomy info for " + aid)

    pid2info_dict = {pid: pid2info_dict.get(pid, {}) for pid in id_list}
    assert len(pid2info_dict) == len(set(id_list))
    access_intermedia(pid2info_dict, suffix=suffix)
    return pid2info_dict


def main(infile, ofile, force=False, redo=False):
    order_id_list, id2annotate = parse_id(infile)
    id2gi = {}
    if isinstance(id2annotate[order_id_list[0]], dict):
        # it is a dict, so it contains other infomation or implemented GI. it may be passed over.
        if 'GI' in id2annotate[order_id_list[0]]:
            print("provided file already contains `GI` column(doesn't check the validation/completeness). Giving `force` param to overwrite/implement it. ")
            if not force:
                id2gi = {k: id2annotate[k]['GI'] for k in order_id_list}
        # todo: re-implemented original infomation into `ofile` from `infile`
    else:
        # no header, just a list of IDs
        pass
    if not id2gi:
        id2gi = pid2GI(order_id_list, redo=redo)

    pid2info_dict = GI2tax(id2gi, redo=redo)
    if not exists(dirname(ofile)) and dirname(ofile):
        os.makedirs(dirname(ofile))

    with open(ofile, 'w') as f1:
        print('#accession ID\tGI\ttaxid\t' + '\t'.join(taxons), file=f1)
        for id, info_dict in pid2info_dict.items():
            GI = info_dict.get('GI', '')
            taxid = info_dict.get('taxid', '')
            print(f"{id}\t{GI}\t{taxid}\t" + '\t'.join([info_dict.get(tax, '') for tax in taxons]), file=f1)
    tqdm.write('finish writing into ' + ofile)


@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id ')
@click.option('-o', 'ofile', help='output file')
@click.option('-f', 'force', help='force overwrite?', default=False, required=False, is_flag=True)
@click.option('-redo', 'redo', help='use cache or not? default is use the cache.', default=False, required=False, is_flag=True)
@click.option('-s', 'start_at', help='start from `protein` or `genome` ID.  etc, protein id maybe like `CBH97221.1`. genome ID should like `GCF_900176205.1` ', default='protein',
              required=False)
def cli(infile, ofile, force, redo,start_at):
    main(infile, ofile, force, redo)


if __name__ == "__main__":
    cli()
