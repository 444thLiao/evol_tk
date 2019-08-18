#!/usr/bin/env python

import os
from subprocess import check_call
from glob import glob
from os.path import join, dirname, abspath
import click

program = "/home-user/software/fs_4.0.1/makeuniformrecfile.pl"


@click.command()
@click.option("-i", "indir", help='input directory which contains `bialle_SNP.ref_*.phas`', default='.', required=False)
@click.option("-log", "log")
def main(indir, log):
    odir = dirname(log)
    os.makedirs(odir, exist_ok=True)
    log_stream = open(log, 'w')
    indir = abspath(indir)
    for ffile in glob(join(indir, "bialle_SNP.ref_*.phase")):
        phasefile = ffile
        recfile = ffile.replace(".phase", ".recfile")
        check_call(f"{program} {phasefile} {recfile}",
                   shell=True,
                   stdout=log_stream,
                   stderr=log_stream)


if __name__ == '__main__':
    main()
