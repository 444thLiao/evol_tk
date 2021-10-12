"""
Advanced script for who wants to modify or manipulate the batch_run script

example use it as batch_mafft

"""

import multiprocessing as mp
import os
from glob import glob
from os.path import *
from subprocess import check_call

import click
from tqdm import tqdm

command_template = "mafft --maxiterate 1000 --genafpair --thread -1 {infile} > {ofile}"


def run(cmd):
    check_call(cmd, shell=True)


def quiet_run(cmd):
    check_call(
        cmd, shell=True, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w")
    )


def main(
    indir,
    odir,
    num_parellel,
    suffix="",
    new_suffix="",
    force=False,
    cmd=command_template,
    test=False,
    get_name=False,
    quiet=False,
):
    suffix = suffix.strip(".")
    new_suffix = new_suffix.strip(".")
    odir = abspath(odir)
    if not exists(odir):
        os.makedirs(odir)
    if suffix:
        suffix = "." + suffix
    file_list = glob(join(indir, f"*{suffix}"))
    file_list = [realpath(abspath(_)) for _ in file_list]
    if not file_list:
        exit(f"empty files, please check your suffix ({indir}/{suffix}) ")
    tqdm.write(
        "start to process %s file with '%s' as suffix" % (len(file_list), suffix)
    )
    params = []
    for infile in tqdm(file_list):
        if new_suffix is not None and suffix is not None:
            if new_suffix == "":
                new_suffix = ""
            else:
                new_suffix = f".{new_suffix.strip('.')}"
                suffix = f".{suffix.strip('.')}"
            ofile = join(odir, basename(infile).replace(suffix, new_suffix))
        else:
            ofile = join(odir, basename(infile))
        if not exists(ofile) or force:
            if get_name:
                name = basename(infile).replace(f"{suffix}", "")
                filled_cmd = cmd.format(infile=infile, ofile=ofile, name=name)
            else:
                filled_cmd = cmd.format(infile=infile, ofile=ofile)
            params.append(filled_cmd)
    if test:
        print(params)
        return
    if quiet:
        with mp.Pool(processes=num_parellel) as tp:
            r = list(tqdm(tp.imap(quiet_run, params), total=len(params)))
    else:
        with mp.Pool(processes=num_parellel) as tp:
            r = list(tqdm(tp.imap(run, params), total=len(params)))


@click.command(
    help="This script accept input directory(-i) which contains files with suffix(-s) and output directory(-o) which will stodge result with its name and new suffix (-ns). It could auto parellel your command into (-np) times. "
)
@click.option("-i", "indir", help="input directory for iterations. ")
@click.option("-o", "odir", help="ouput directory for stodge the output files")
@click.option(
    "-s",
    "suffix",
    default=None,
    help="suffix of input files needed to be iterated within the indir,default is empty",
)
@click.option(
    "-ns",
    "new_suffix",
    default=None,
    help="new suffix of output files, default is empty",
)
@click.option(
    "-np",
    "num_parellel",
    default=10,
    help="num of processes could be parellel.. default is 10",
)
@click.option(
    "-f",
    "force",
    default=False,
    required=False,
    is_flag=True,
    help="overwrite the output files or not.",
)
@click.option(
    "-quiet",
    "quiet",
    default=False,
    required=False,
    is_flag=True,
    help="overwrite the output files or not.",
)
@click.option("-t", "test", help="test?", default=False, required=False, is_flag=True)
@click.option(
    "-get_name",
    "get_name",
    default=False,
    required=False,
    is_flag=True,
    help="get the basename of the input file as a new format string to the command. Use {name} in the command to use it. ",
)
@click.option(
    "-cmd",
    "cmd",
    help="it shoulw accept a command with {} as indicator of string format. e.g. mafft --maxiterate 1000 --genafpair --thread -1 {infile} > {ofile}, the suffix of original file and new file could be ignore and it will auto added by the script. The suffix should be assigned at parameter `ns` or `s`. Default both are empty. If you want to add more flexible parameters, it should modify this script directly. Beside that, the `get_name` parameter could help you to extract the basename of the input file and pass it to the `cmd`. ",
)
def cli(
    indir, odir, suffix, new_suffix, force, test, num_parellel, cmd, get_name, quiet
):
    main(
        indir=indir,
        odir=odir,
        num_parellel=num_parellel,
        suffix=suffix,
        new_suffix=new_suffix,
        force=force,
        cmd=cmd,
        test=test,
        get_name=get_name,
        quiet=quiet,
    )


if __name__ == "__main__":
    cli()
