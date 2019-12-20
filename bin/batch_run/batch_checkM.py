"""
For parellel checkm, after that you could merge it by yourself
"""

import multiprocessing as mp
import os
from glob import glob
from os.path import *
# sys.path.insert(0,dirname(dirname(dirname(dirname(__file__)))))
from subprocess import check_call

import click
from tqdm import tqdm

command_template = 'checkm taxonomy_wf -t 10 {extra_option} -x {infile} {tax} {indir} {odir}'


def run(args):
    unit_run(*args)


def unit_run(infile, indir, tax, odir,extra_option):
    check_call(command_template.format(infile=infile,
                                       indir=indir,
                                       tax=tax,
                                       odir=odir,
                                       extra_option=extra_option),
               shell=True,
               stdout=open('/dev/null', 'w'))


def main(indir, odir, tax, use_fa,num_parellel, suffix='', force=False):
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
                           tax,
                           new_odir,
                           '' if use_fa else '-g' ))
    with mp.Pool(processes=num_parellel) as tp:
        r = list(tqdm(tp.imap(run, params), total=len(params)))


@click.command()
@click.option('-i', 'indir')
@click.option('-o', 'odir')
@click.option('-s', 'suffix', default='')
@click.option('-np', 'num_parellel', default=10)
@click.option('-t', 'tax', default='domain Bacteria')
@click.option('-f', 'force', help='overwrite?', default=False, required=False, is_flag=True)
@click.option('-use_fa', 'use_fa', help='use nucleotide sequence or not.default using annotated proteins?', default=False, required=False, is_flag=True)
def cli(indir, odir, tax, suffix, force, num_parellel,use_fa):
    main(indir=indir,
         odir=odir,
         tax=tax,
         num_parellel=num_parellel,
         suffix=suffix,
         force=force,
         use_fa=use_fa)


if __name__ == "__main__":
    cli()
