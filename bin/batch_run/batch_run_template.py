import multiprocessing as mp
import os
from glob import glob
from os.path import *
# sys.path.insert(0,dirname(dirname(dirname(dirname(__file__)))))
from subprocess import check_call

import click
from tqdm import tqdm
from dating_workflow.step_script import get_files

command_template = 'mafft --maxiterate 1000 --genafpair --thread -1 {infile} > {ofile}'


def run(args):
    unit_run(*args)


def unit_run(infile, ofile):
    check_call(command_template.format(infile=infile,
                                       ofile=ofile),
               shell=True)


def main(indir, odir, num_parellel, suffix='', new_suffix='', force=False):
    suffix = '.'+suffix.strip('.')
    new_suffix = '.'+new_suffix.strip('.')
    if not exists(odir):
        os.makedirs(odir)
    file_list = get_files(indir,suffix)
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
            params.append((in_file, ofile))
    with mp.Pool(processes=num_parellel) as tp:
        r = list(tqdm(tp.imap(run, params), total=len(params)))


@click.command()
@click.option('-i', 'indir')
@click.option('-o', 'odir')
@click.option('-s', 'suffix', default='')
@click.option('-ns', 'new_suffix', default='')
@click.option('-np', 'num_parellel', default=10)
@click.option('-f', 'force', help='overwrite?', default=False, required=False, is_flag=True)
def cli(indir, odir, suffix, new_suffix, force, num_parellel):
    main(indir=indir,
         odir=odir,
         num_parellel=num_parellel,
         suffix=suffix,
         new_suffix=new_suffix,
         force=force)


if __name__ == "__main__":
    cli()
