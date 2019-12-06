"""
For parellel checkm, after that you could merge it by yourself
"""

import sys
from os.path import *
import os
# sys.path.insert(0,dirname(dirname(dirname(dirname(__file__)))))
from subprocess import check_call
import click
from glob import glob
from tqdm import tqdm
import multiprocessing as mp

command_template = 'checkm taxonomy_wf -t 10 -x {infile} domain Bacteria {indir} {odir}'


def run(args):
    unit_run(*args)


def unit_run(infile,indir, odir):
    check_call(command_template.format(infile=infile,
                                       indir=indir,
                                       odir=odir),
               shell=True,
               stdout=open('/dev/null','w'))


def main(indir, odir, num_parellel, suffix='', force=False):
    suffix = suffix.strip('.')
    if not exists(odir):
        os.makedirs(odir)
    if suffix:
        suffix = '.' + suffix
    file_list = glob(join(indir, f'*{suffix}'))
    tqdm.write("start to process %s file with '%s' as suffix" % (len(file_list), suffix))
    params = []
    for in_file in file_list:
        new_odir = join(odir,
                        basename(in_file).replace(f'{suffix}',
                                                       ''))
        if not exists(new_odir) or force:
            params.append((basename(in_file),
                           indir, 
                           new_odir))
    with mp.Pool(processes=num_parellel) as tp:
        r = list(tqdm(tp.imap(run,params),total=len(params)))


@click.command()
@click.option('-i', 'indir')
@click.option('-o', 'odir')
@click.option('-s', 'suffix', default='')
@click.option('-np', 'num_parellel', default=10)
@click.option('-f', 'force', help='overwrite?', default=False, required=False, is_flag=True)
def cli(indir, odir, suffix, new_suffix, force, num_parellel):
    main(indir=indir,
         odir=odir,
         num_parellel=num_parellel,
         suffix=suffix,
         force=force)


if __name__ == "__main__":
    cli()
