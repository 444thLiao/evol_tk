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
    convertor.get_taxon()

    pid2info_dict = defaultdict(dict)
    for oid in convertor.origin_ids:
        pid2info_dict[oid]['GI'] = convertor.GI[oid]
        pid2info_dict[oid]['taxid'] = convertor.tids[oid]
        pid2info_dict[oid]['accession'] = oid
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

    # if start_at == 'protein':
    #     if isinstance(id2annotate[order_id_list[0]], dict):
    #         # it is a dict, so it contains other information or implemented GI. it may be passed over.
    #         if 'GI' in id2annotate[order_id_list[0]]:
    #             print("provided file already contains `GI` column(doesn't check the validation/completeness). Giving `force` param to overwrite/implement it. ")
    #             if not force:
    #                 id2gi = {k: id2annotate[k]['GI'] for k in order_id_list}
    #
    #     else:
    #         # no header, just a list of IDs
    #         pass
    # elif start_at == 'genome':
    #     exit("Use genome2tax instead")
    # else:
    #     raise SyntaxError('wrong input of start_at')


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
