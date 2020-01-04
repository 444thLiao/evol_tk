import multiprocessing as mp
import os
from glob import glob
from os.path import *
# sys.path.insert(0,dirname(dirname(dirname(dirname(__file__)))))
from subprocess import check_call

import click
from tqdm import tqdm

command_template = "trimal -in {in_file} -out {o_file} -automated1 -resoverlap {resoverlap} -seqoverlap {seqoverlap}"


def run(args):
    unit_run(*args)


def unit_run(in_file, o_file,resoverlap,seqoverlap):
    check_call(command_template.format(in_file=in_file,
                                       o_file=o_file,
                                       resoverlap=resoverlap,
                                       seqoverlap=seqoverlap),
               shell=1)


def main(in_dir, odir, num_parellel, suffix='', new_suffix='',resoverlap=0.55,seqoverlap=60,**kwarg):
    suffix = suffix.strip('.')
    new_suffix = new_suffix.strip('.')
    if not exists(odir):
        os.makedirs(odir)
    if suffix:
        suffix = '.' + suffix
    file_list = glob(join(in_dir, f'*{suffix}'))
    tqdm.write("start to process %s file with '%s' as suffix" % (len(file_list), suffix))
    params = []
    for in_file in file_list:
        if new_suffix and suffix:
            ofile = join(odir,
                         basename(in_file).replace(suffix,
                                                   '.' + new_suffix))
        else:
            ofile = join(odir,
                         basename(in_file))
        params.append((in_file, ofile,resoverlap,seqoverlap))
    with mp.Pool(processes=num_parellel) as tp:
        r = list(tqdm(tp.imap(run, params), total=len(params)))


@click.command()
@click.option('-i', 'indir')
@click.option('-o', 'odir')
@click.option('-s', 'suffix', default='aln')
@click.option('-ns', 'new_suffix', default='trimal')
@click.option('-np', 'num_parellel', default=10)
@click.option('-ro','resoverlap',default=0.55)
@click.option('-so','seqoverlap',default=60)
def cli(indir, odir, num_parellel, suffix, new_suffix,resoverlap,seqoverlap):
    main(indir, odir, num_parellel, suffix, new_suffix,resoverlap,seqoverlap)


if __name__ == "__main__":
    cli()
