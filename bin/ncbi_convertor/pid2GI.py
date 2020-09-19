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

from bin.ncbi_convertor import edl, access_intermedia, parse_id,NCBI_convertor,process_path

def pid2GI(id_list, redo=False):
    convertor = NCBI_convertor(id_list,"protein")
    suffix = 'protein2GI'
    convertor.check_cache(suffix=suffix,redo=redo)

    convertor.get_GI()
    id2gi = convertor.GI
    return id2gi


def main(infile, ofile, force=False, redo=False):
    ofile = process_path(ofile)
    order_id_list, id2annotate = parse_id(infile)
    id2gi = pid2GI(order_id_list, redo=redo)

    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))

    with open(ofile, 'w') as f1:
        print('#accession ID\tGI', file=f1)
        for id, GI in id2gi.items():
            print(f'{id}\t{GI}', file=f1)
    tqdm.write('finish writing into ' + ofile)

    # too verbose
    # if isinstance(id2annotate[order_id_list[0]], dict):
    #     # it is a dict, so it contains other infomation or implemented GI. it may be passed over.
    #     if 'GI' in id2annotate[order_id_list[0]]:
    #         print("provided file already contains `GI` column(doesn't check the validation/completeness). Giving `force` param to overwrite/implement it. ")
    #         if not force:
    #             id2gi = {k: id2annotate[k]['GI'] for k in order_id_list}
    # else:
    #     # no header, just a list of IDs
    #     pass



@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id ')
@click.option('-o', 'ofile', help='output file')
@click.option('-f', 'force', help='force overwrite?', default=False, required=False, is_flag=True)
@click.option('-redo', 'redo', help='use cache or not? default is use the cache.', default=False, required=False, is_flag=True)
def cli(infile, ofile, force, redo):
    main(infile, ofile, force, redo)


if __name__ == "__main__":
    cli()
