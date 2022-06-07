"""
For easily perform hmmscan
"""

import multiprocessing as mp
import os
from glob import glob
from os.path import *
from subprocess import check_call

import click
from Bio import SeqIO
from tqdm import tqdm
from dating_workflow.step_script import get_files
default_db = '/home-db/pub/protein_db/TIGRFAM/14.0_release/TIGRFAM.HMM'
command_template = '/home-user/thliao/anaconda3/bin/hmmscan --tblout {o_file} --acc --noali --notextw --cpu 20 {db} {in_file} > /dev/null 2>&1'


def run(args):
    unit_run(*args)


def unit_run(in_file, o_file, db):
    check_call(command_template.format(in_file=in_file,
                                       o_file=o_file,
                                       db=db),
               shell=True,
               stdout=open('/dev/null', 'w'),
               stderr=open('/dev/null', 'w'))

from api_tools.tk import convert_genome_ID_rev,convert_genome_ID


def main(in_dir, odir, num_parellel, suffix='', gids=None, force=False, db=default_db, **kwarg):
    suffix = '.' + suffix.strip('.')
    new_suffix = '.tab'
    if not exists(odir):
        os.makedirs(odir)
    file_list = get_files(in_dir,suffix)
    if gids is not None:
        gids = set(gids)
        os.makedirs(join(odir, 'tmp'), exist_ok=1)
        new_file_list = []
        tqdm.write('iterating files to collect with giving genome ids')
        for f in tqdm(file_list):
            records = SeqIO.parse(f, format='fasta')
            records = [_
                       for _ in records
                       if _.id in gids]
            if not records:
                records = SeqIO.parse(f, format='fasta')
                records = [_
                           for _ in records
                           if convert_genome_ID_rev(_.id.split('_')[0]) in gids]
            n_f = join(odir, 'tmp', basename(f))
            if not records:
                raise Exception('error not records')
            with open(n_f, 'w') as f1:
                SeqIO.write(records, f1, format='fasta-2line')
            new_file_list.append(n_f)
        file_list = new_file_list[::]
    tqdm.write("start to process %s file with '%s' as suffix" % (len(file_list), suffix))
    params = []
    for in_file in tqdm(file_list):
        if new_suffix and suffix:
            ofile = join(odir,
                         basename(in_file).replace(suffix,
                                                   new_suffix))
        else:
            ofile = join(odir,
                         basename(in_file))
        if not exists(ofile) or force:
            params.append((in_file, ofile, db))
    with mp.Pool(processes=num_parellel) as tp:
        r = list(tqdm(tp.imap(run, params), total=len(params)))


@click.command()
@click.option('-i', 'indir')
@click.option('-o', 'odir')
@click.option('-s', 'suffix', help='suffix for files', default='faa')
@click.option('-np', 'num_parellel', default=10)
@click.option("-gl", "genome_list", default=None,
              help="It could be None. If you provided, you could use it to subset the files by indicate names.")
@click.option('-db', 'db', default=default_db)
@click.option('-f', 'force', help='overwrite?', default=False, required=False, is_flag=True)
def cli(indir, odir, num_parellel, suffix, genome_list, force, db):
    if genome_list is None:
        gids = None
    else:
        gids = open(genome_list).read().split('\n')
        gids = set([_ for _ in gids if _])
    main(indir, odir, num_parellel, suffix, gids=gids, force=force, db=db)


if __name__ == "__main__":
    cli()
