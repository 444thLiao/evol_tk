"""
This script is mainly for implementing genome ID and its assembly ID with its protein accession id.
It will generate a tab separated table.
Because each protein ID would corresponding a lots of identical protein groups (which this script mainly used database)
it will generate more line than input ID. but the set of/deduplicated column 1 must match the input IDs. If not ,it will leave it blank row.

For filtration amplicons, it could simply find proteins which only contain NONE assembly ID columns. (it simple!)
header is ['accession ID', 'GI', 'assembly_ID', 'nuccore ID', 'start', 'end', 'strand']
or look like accession ID\tGI\tassembly_ID\tnuccore ID\tstart\tend\tstrand
"""
import os
from os.path import exists, dirname

import click
from tqdm import tqdm

from bin.ncbi_convertor import parse_id, NCBI_convertor, process_path

def main(infile, ofile, force=False):
    ofile = process_path(ofile)
    order_id_list, id2annotate = parse_id(infile)
    convertor = NCBI_convertor(order_id_list,db='protein')
    pid2assembly_dict = convertor.pid2assembly()

    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))

    if exists(ofile) and not force:
        tqdm.write("detect existing " + ofile + ' no force param input, so it quit instead of writing.')
        return

    with open(ofile, 'w') as f1:
        print('#accession ID\tGI\tassembly_ID\tnuccore ID\tstart\tend\tstrand', file=f1)
        for pid in convertor.origin_ids:
            if convertor.GI is not None:
                GI = convertor.GI[pid]
            else:
                GI = 'NA'
            _dict = pid2assembly_dict[pid]
            for assembly_id,info in zip(_dict['assembly'],_dict['nuc_info']):
                print(f'{pid}\t{GI}\t{assembly_id}\t' + '\t'.join(map(str, info)), file=f1)
    tqdm.write('finish writing into ' + ofile)


@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id ')
@click.option('-o', 'ofile', help='output file')
@click.option('-f', 'force', help='force overwrite?', default=False, required=False, is_flag=True)
def cli(infile, ofile, force):
    main(infile, ofile, force)


if __name__ == "__main__":
    cli()
