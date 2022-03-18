# -*- coding: utf-8 -*-
"""
This script is mainly for implementing taxonomy information with its protein accession id.
It will generate a tab separated table.
"""
import os
from collections import defaultdict
from os.path import exists, dirname

import click
from tqdm import tqdm

from bin.ncbi_convertor import access_intermedia, parse_id, taxons, NCBI_convertor, process_path


def id2tax(id_list, redo=False, db="protein"):
    convertor = NCBI_convertor(id_list, db)
    suffix = 'protein2GI'
    convertor.check_cache(suffix=suffix, redo=redo)
    tids = convertor.get_taxon('update')

    pid2info_dict = defaultdict(dict)
    for oid in convertor.origin_ids:
        if convertor.GI is not None:
            # NCBI has removed the GI gradually. We could use the id directly instead of transfering it into GI.
            pid2info_dict[oid]['GI'] = convertor.GI[oid]
        pid2info_dict[oid]['taxid'] = tids.get(oid,'NA')
        pid2info_dict[oid]['accession'] = oid
    id2taxon = convertor.get_taxons_from_tid(tids)
    for pid,_d in pid2info_dict.items():
        _d.update(id2taxon.get(pid,{}))
    for id in id_list:
        if id not in pid2info_dict:
            pid2info_dict[id] = {}
    access_intermedia(pid2info_dict, suffix=suffix)
    return pid2info_dict


def main(infile, ofile, db='protein', redo=False):
    ofile = process_path(ofile)
    order_id_list, id2annotate = parse_id(infile)
    pid2info_dict = id2tax(order_id_list, redo=redo, db=db)
    if not exists(dirname(ofile)):
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
@click.option('-db', 'source_db', help='start from `protein` or `assembly` ID.  etc, protein id maybe like `CBH97221.1`. assembly ID should like `GCF_900176205.1` ',
              default='protein',
              required=False)
def cli(infile, ofile, source_db):
    source_db = source_db.strip().lower()
    if source_db not in ['assembly', 'protein']:
        raise IOError('Unexpected params of start_at. giving `%s`?? ' % source_db)
    main(infile, ofile, db=source_db)


if __name__ == "__main__":
    cli()
