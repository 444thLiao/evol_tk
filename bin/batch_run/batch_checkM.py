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

from dating_workflow.step_script import get_files

command_template = '`which checkm` taxonomy_wf -t 10 {extra_option} -x {infile} {tax} {indir} {odir}'


def run(args):
    unit_run(*args)


def unit_run(infile, indir, tax, odir, extra_option):
    check_call(command_template.format(infile=infile,
                                       indir=indir,
                                       tax=tax,
                                       odir=odir,
                                       extra_option=extra_option),
               shell=True,
               stdout=open('/dev/null', 'w'))


def main(indir, odir, tax, use_fa, num_parellel, suffix='', force=False, genome_list=None):
    suffix = '.'+suffix.strip('.') # format it
    if not exists(odir):
        os.makedirs(odir)
    file_list = get_files(indir,suffix)
    
    subset_names = []
    if genome_list is not None:
        subset_names = [_ for _ in open(genome_list).read().split('\n') if _]

    # suffix here already have . as prefix.
    if subset_names:
        file_list = [_
                     for _ in file_list
                     if basename(_).replace(suffix, '') in subset_names]
    else:
        pass
    tqdm.write("start to process %s file with '%s' as suffix" % (len(file_list),
                                                                 suffix))
    params = []

    for in_file in file_list:
        gname = basename(in_file).replace(suffix, '')
        new_odir = join(odir,
                        gname)
        if not exists(new_odir) or force:
            params.append((basename(in_file),
                           indir,
                           tax,
                           new_odir,
                           '' if use_fa else '-g'))
    with mp.Pool(processes=num_parellel) as tp:
        r = list(tqdm(tp.imap(run, params), total=len(params)))


@click.command(help="batch run pipelines 'taxonomy_wf' of checkm.")
@click.option('-i', 'indir', help='input directory')
@click.option('-o', 'odir', help='output directory')
@click.option('-s', 'suffix', help='suffix of used files', default='')
@click.option('-np', 'num_parellel', default=10)
@click.option('-t', 'tax', default='domain Bacteria', help='designated taxonomic level you want to use.')
@click.option('-f', 'force', help='overwrite?', default=False, required=False, is_flag=True)
@click.option('-use_fa', 'use_fa', help='use nucleotide sequence or not. default using annotated proteins?', default=False, required=False, is_flag=True)
@click.option("-gl", "genome_list", default=None,
              help="It could be None. If you provided, you could use it to subset the used sequences by its name.")
def cli(indir, odir, tax, suffix, force, num_parellel, use_fa, genome_list):
    main(indir=indir,
         odir=odir,
         tax=tax,
         num_parellel=num_parellel,
         suffix=suffix,
         force=force,
         use_fa=use_fa,
         genome_list=genome_list)


if __name__ == "__main__":
    cli()
