import multiprocessing as mp
import os
from glob import glob
from os.path import *
# sys.path.insert(0,dirname(dirname(dirname(dirname(__file__)))))
from subprocess import check_call

import click
from tqdm import tqdm

command_template = 'iqtree -nt 20 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre {ofile} -s {infile}'
command_template2 = '/home-user/thliao/anaconda3/bin/FastTreeMP {infile} > {ofile}'

used_command = {'iqtree': command_template,
                'fasttree': command_template2}


def run(args):
    unit_run(*args)


def unit_run(infile, ofile, cmd):
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    check_call(cmd.format(infile=infile,
                          ofile=ofile),
               shell=True,
               stdout=open('/dev/null', 'w'),
               stderr=open('/dev/null', 'w'))


def main(indir, odir, num_parellel, suffix='', new_suffix='', force=False, software='iqtree'):
    suffix = suffix.strip('.')
    new_suffix = new_suffix.strip('.')
    if software not in used_command:
        raise IOError
    cmd = used_command[software]
    if software == 'fasttree' and new_suffix == 'iqtree':
        new_suffix = 'newick'
    if not exists(odir):
        os.makedirs(odir)
    if suffix:
        suffix = '.' + suffix
    file_list = glob(join(indir, f'*{suffix}'))
    tqdm.write("start to process %s file with '%s' as suffix" % (len(file_list), suffix))
    params = []
    for in_file in tqdm(file_list):
        name = basename(in_file).replace(suffix, '')
        if new_suffix and suffix:
            ofile = join(odir,
                         name,
                         basename(in_file).replace(suffix,
                                                   '.' + new_suffix))
        else:
            ofile = join(odir,
                         name,
                         basename(in_file))
        if not exists(ofile) or force:
            params.append((in_file, ofile, cmd))
    with mp.Pool(processes=num_parellel) as tp:
        r = list(tqdm(tp.imap(run, params), total=len(params)))


@click.command()
@click.option('-i', 'indir')
@click.option('-o', 'odir')
@click.option('-s', 'suffix', default='aln')
@click.option('-ns', 'new_suffix', default='iqtree')
@click.option('-np', 'num_parellel', default=10)
@click.option('-f', 'force', help='overwrite?', default=False, required=False, is_flag=True)
@click.option('-use', 'software', help='which software used to build tree?', default='iqtree')
def cli(indir, odir, suffix, new_suffix, force, num_parellel, software):
    main(indir=indir,
         odir=odir,
         num_parellel=num_parellel,
         suffix=suffix,
         new_suffix=new_suffix,
         force=force,
         software=software)


if __name__ == "__main__":
    cli()
