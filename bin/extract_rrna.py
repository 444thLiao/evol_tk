"""
From prokka output to extract rrna (mainly for 16S and 23S)
it follow the structure of traditional prokka output
"""

from Bio import SeqIO
from glob import glob
from os.path import *
from tqdm import tqdm
from collections import defaultdict
import click
import os


def extract_16s(gbk_files):
    name2rrna = defaultdict(list)
    for f in tqdm(gbk_files):
        name = basename(dirname(f))
        if not exists(f):
            print('no sequence in ', name)
            continue
        f = f.replace('.gbk', '.ffn')
        records = SeqIO.parse(f, format='fasta')
        records = list(records)
        name = basename(dirname(f))
        rrna_16S = {str(_.seq): _
                    for _ in records
                    if '16S ribosomal RNA' in _.description}

        rrna_23S = {str(_.seq): _
                    for _ in records
                    if '23S ribosomal RNA' in _.description}

        name2rrna[name] = (list(rrna_16S.values()),
                           list(rrna_23S.values()))
    return name2rrna


def get_stats(name2rrna):
    no16s_count = 0
    no23s_count = 0
    multiple_16S_count = 0
    multiple_23S_count = 0
    for name, (rrna_16S, rrna_23S) in name2rrna.items():
        if not rrna_16S:
            no16s_count += 1
        if not rrna_23S:
            no23s_count += 1
        if len(rrna_16S) > 1:
            multiple_16S_count += 1
        if len(rrna_23S) > 1:
            multiple_23S_count += 1
    print(f"{no16s_count} genomes doesn't contain 16S ")
    print(f"{no23s_count} genomes doesn't contain 23S ")
    print(f"{multiple_16S_count} genomes contains multiple 16S ")
    print(f"{multiple_23S_count} genomes contains multiple 23S ")
    return


def process_path(path):
    if not '/' in path:
        path = './' + path
    if path.startswith('~'):
        path = expanduser(path)
    if path.startswith('.'):
        path = abspath(path)
    return path


@click.command()
@click.option('-i', 'indir', help='directory which is prokka_o. or some special soft link stodge dir')
@click.option('-o', 'odir', help='normal, it will generate 16S and 23S files separately. ')
@click.option('-ps', 'preset', help='', default='prokka')
def main(indir, odir, preset):
    indir = process_path(indir)
    odir = process_path(odir)
    if not exists(odir):
        os.makedirs(odir)

    gbk_files = []
    if preset == 'prokka':
        gbk_files = glob(join(indir, '*', '*.gbk'))
    elif preset == 'dating':
        pfiles = glob(join(indir, '*.faa'))
        all_ids = set([basename(_).replace('.faa', '') for _ in pfiles])
        pfiles = [realpath(_) for _ in pfiles]
        target_d = set([dirname(dirname(_)) for _ in pfiles])
        for _ in target_d:
            gbk_files += list(glob(join(_, 'tmp', '*', '*.gbk')))
        gbk_files = [_ for _ in gbk_files
                     if basename(dirname(_)) in all_ids]
    name2rrna = extract_16s(gbk_files)
    get_stats(name2rrna)
    records_16S = [_
                   for rrna_16S, rrna_23S in name2rrna.values()
                   for _ in rrna_16S]
    records_23S = [_
                   for rrna_16S, rrna_23S in name2rrna.values()
                   for _ in rrna_23S]
    with open(join(odir, '16S.fasta'), 'w') as f1:
        SeqIO.write(records_16S, f1, format='fasta-2line')
    with open(join(odir, '23S.fasta'), 'w') as f1:
        SeqIO.write(records_23S, f1, format='fasta-2line')


if __name__ == '__main__':
    import sys

    #assert len(sys.argv) == 3
    main()
    # extract_rrna.py -i ./pipelines_o/prokka_o -o /home-user/thliao/data/jjtao_20191113/
    # if using dating
    # extract_rrna.py -i ~/data/nitrification_for/dating_for/raw_genome_proteins -o ~/data/nitrification_for/dating_for/raw_genome_proteins -ps dating
