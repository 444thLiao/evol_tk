import multiprocessing as mp
import os
from glob import glob
from os.path import *
# sys.path.insert(0,dirname(dirname(dirname(dirname(__file__)))))
from subprocess import check_call

import click
from tqdm import tqdm
from dating_workflow.step_script import get_files
command_template = "trimal -in {in_file} -out {o_file} -automated1 -resoverlap {resoverlap} -seqoverlap {seqoverlap}"


def run(args):
    unit_run(*args)


def unit_run(in_file, o_file, resoverlap, seqoverlap):
    check_call(command_template.format(in_file=in_file,
                                       o_file=o_file,
                                       resoverlap=resoverlap,
                                       seqoverlap=seqoverlap),
               shell=1)


def main(in_dir, odir, num_parellel, suffix='', new_suffix='', resoverlap=0.55, seqoverlap=60, **kwarg):
    suffix = '.'+suffix.strip('.')
    new_suffix = '.'+new_suffix.strip('.')
    if not exists(odir):
        os.makedirs(odir)
    file_list = get_files(in_dir,suffix)
    tqdm.write("start to process %s file with '%s' as suffix" % (len(file_list), suffix))
    params = []
    for in_file in file_list:
        if new_suffix and suffix:
            ofile = join(odir,
                         basename(in_file).replace(suffix,
                                                   new_suffix))
        else:
            ofile = join(odir,
                         basename(in_file))
        params.append((in_file, ofile, resoverlap, seqoverlap))
    with mp.Pool(processes=num_parellel) as tp:
        r = list(tqdm(tp.imap(run, params), total=len(params)))


@click.command(help="")
@click.option('-i', 'indir',help="input directory which stodge files with suffix. Normally, it should contain aln files directly.  ")
@click.option('-o', 'odir',default=None,help="output directory which stodge the ouput aln files. default is the identifical to the input dir. Inplaced output. ")
@click.option('-s', 'suffix', default='aln',help="suffix")
@click.option('-ns', 'new_suffix', default='trimal',help="new suffix")
@click.option('-np', 'num_parellel', default=10,help="number of parellel processes you want to run. ")
@click.option('-ro', 'resoverlap', default=0.55,help="Minimum overlap of a positions with other positions in the column to be considered a good position. Range: [0 - 1]. See the userguide of trimal ")
@click.option('-so', 'seqoverlap', default=60,help="Minimum percentage of good positions that a sequence must have in order to be conserved. Range: [0 - 100]. See the userguide of trimal")
def cli(indir, odir, num_parellel, suffix, new_suffix, resoverlap, seqoverlap):
    if odir is None:
        odir = indir
    main(indir, odir, num_parellel, suffix, new_suffix, resoverlap, seqoverlap)


if __name__ == "__main__":
    cli()
