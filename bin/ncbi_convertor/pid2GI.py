"""
This script is the first kit of ncbi convertor. It is also a simple example for this convertor.
Other convertor would follow the structure of this script.
"""
import io
import os
from os.path import exists, dirname

import click
from Bio import Entrez
from tqdm import tqdm

from bin.ncbi_convert import edl, access_intermedia, parse_id


def pid2GI(id_list, redo=False):
    suffix = 'pid2gi'
    _cache = access_intermedia(id_list, suffix=suffix, redo=redo)
    if _cache is not None:
        id2gi = _cache
    else:
        tqdm.write('from protein Accession ID to GI')
        results, failed = edl.esearch(db='protein',
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

    if not exists(dirname(ofile)) and dirname(ofile):
        os.makedirs(dirname(ofile))

    with open(ofile, 'w') as f1:
        print('#accession ID\tGI', file=f1)
        for id, GI in id2gi.items():
            print(f'{id}\t{GI}', file=f1)
    tqdm.write('finish writing into ' + ofile)


@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id ')
@click.option('-o', 'ofile', help='output file')
@click.option('-f', 'force', help='force overwrite?', default=False, required=False, is_flag=True)
@click.option('-redo', 'redo', help='use cache or not? default is use the cache.', default=False, required=False, is_flag=True)
def cli(infile, ofile, force, redo):
    main(infile, ofile, force, redo)


if __name__ == "__main__":
    cli()
