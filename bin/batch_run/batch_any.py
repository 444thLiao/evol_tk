"""
Advanced script for who wants to modify or manipulate the batch_run script
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

command_template = 'mafft --maxiterate 1000 --genafpair --thread -1 {infile} > {ofile}'


def run(cmd):
    check_call(cmd,
               shell=True)


def unit_run(infile, ofile):
    check_call(command_template.format(infile=infile,
                                       ofile=ofile),
               shell=True)


def main(indir, odir, num_parellel, suffix='', new_suffix='', force=False,cmd=command_template):
    suffix = suffix.strip('.')
    new_suffix = new_suffix.strip('.')
    odir = abspath(odir)
    if not exists(odir):
        os.makedirs(odir)
    if suffix:
        suffix = '.' + suffix
    file_list = glob(join(indir, f'*{suffix}'))
    if not file_list:
        exit(f"empty files, please check your suffix ({indir}/{suffix}) ")
    tqdm.write("start to process %s file with '%s' as suffix" % (len(file_list), suffix))
    params = []
    for infile in tqdm(file_list):
        if new_suffix and suffix:
            ofile = join(odir,
                         basename(infile).replace(suffix,
                                                   '.' + new_suffix))
        else:
            ofile = join(odir,
                         basename(infile))
        if not exists(ofile) or force:
            cmd = cmd.format(infile=infile,
                                       ofile=ofile)
            params.append(cmd)
    with mp.Pool(processes=num_parellel) as tp:
        r = list(tqdm(tp.imap(run,params),total=len(params)))


@click.command(help="This script accept input directory(-i) which contains files with suffix(-s) and output directory(-o) which will stodge result with its name and new suffix (-ns). It could auto parellel your command into (-np) times. ")
@click.option('-i', 'indir')
@click.option('-o', 'odir')
@click.option('-s', 'suffix', default='')
@click.option('-ns', 'new_suffix', default='')
@click.option('-np', 'num_parellel', default=10)
@click.option('-f', 'force', help='overwrite?', default=False, required=False, is_flag=True)
@click.option('-cmd',"cmd",help="it shoulw accept a command with {} as indicator of string format. e.g. mafft --maxiterate 1000 --genafpair --thread -1 {infile} > {ofile}, the suffix of original file and new file could be ignore. The suffix should be assigned at parameter `ns` or `s`. now default is empty. If you want to add more flexible parameters, it should modify this script directly. ")
def cli(indir, odir, suffix, new_suffix, force, num_parellel,cmd):
    main(indir=indir,
         odir=odir,
         num_parellel=num_parellel,
         suffix=suffix,
         new_suffix=new_suffix,
         force=force,
         cmd=cmd)


if __name__ == "__main__":
    cli()
