"""
from fasta to get id,name into a list
"""
import sys
from os.path import dirname

sys.path.insert(0, dirname(__file__))
import click
# required format for binary python
import os
from os.path import *
from Bio import SeqIO
from tqdm import tqdm


def main(infa, ofile, from_nr):
    records = SeqIO.parse(infa, format='fasta')
    if '/' not in ofile:
        ofile = './' + ofile
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    with open(ofile, 'w') as f1:

        for record in tqdm(records):
            if from_nr:
                ids = record.description.split('\x01')
                for id in ids:
                    id = id.split(' ')[0]
                    f1.write('\t'.join([id, record.name]) + '\n')
            else:
                f1.write('\t'.join([record.id, record.name]) + '\n')


@click.command()
@click.option('-i', 'infa')
@click.option('-o', 'ofile')
@click.option('-nr', 'from_nr', help='come from nr or not',
              default=False, required=False, is_flag=True)
def cli(infa, ofile, from_nr):
    main(infa, ofile, from_nr)


if __name__ == "__main__":
    cli()
