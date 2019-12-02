
from ete3 import Tree
import click
from glob import glob
import multiprocessing as mp
from subprocess import check_call, check_output
import subprocess
import os
from os.path import *
from tqdm import tqdm
indir = '/home-user/thliao/data/nitrification_for/dating_for/mcmc_for/proteins'

indir = '/home-user/thliao/data/nitrification_for/dating_for/mcmc_for/test_no_filled_gaps'
in_phyfile = './concat_aln.phy'
in_treefile = './243g_120gene.3calibrations.newick'
aaRatefile = '../lg.dat'

def modify(file, **kwargs):
    text = open(file).read()
    text = text.split('\n')
    new_text = []
    for row in text:
        key = row.split('=')[0].strip()
        if key in kwargs:
            new_text.append(f"{key} = {kwargs[key]}")
        else:
            new_text.append(row)
    return '\n'.join(new_text)


def run(args):
    if isinstance(args, str):
        cmd = args
    else:
        cmd, log = args
    try:
        subprocess.check_output(cmd, shell=1)
    except subprocess.CalledProcessError as e:
        print('error')
        error = e.output
        return error


def generate_tmp(template_ctl, ndata, odir):
    # template_ctl_01 = './01_mcmctree.ctl'
    new_01_ctl = join(odir, '01_mcmctree_modify.ctl')
    params = {'ndata': ndata,
              'seqfile': abspath(in_phyfile),
              'treefile': abspath(in_treefile),
              'outfile': './01_out'}
    text = modify(template_ctl, **params)
    with open(new_01_ctl, 'w') as f1:
        f1.write(text)
    run(f"export PATH=''; cd {odir}; /home-user/thliao/software/paml4.9j/bin/mcmctree 01_mcmctree_modify.ctl")


def collecting_tmp(tmp_indir, odir):
    if not exists(odir):
        os.makedirs(odir)
    ctls = glob(join(tmp_indir, 'tmp0*.ctl'))
    for ctl in ctls:
        name = basename(ctl).replace('.ctl', '')
        os.makedirs(join(odir, './{name}'), exist_ok=1)
        os.system(f'mv {name}.* {name}/')


def run_each_tmp(tmp_indir, odir, aaRatefile=aaRatefile):
    if not exists(odir):
        os.makedirs(odir)
    params = []
    ctls = glob(join(tmp_indir, '*', 'tmp0*.ctl'))
    ctls = [_ for _ in ctls if 'modify' not in _]  # remove modify
    for ctl in ctls:
        new_text = modify(ctl,
                          **{'model': 2,
                             'aaRatefile': abspath(aaRatefile),
                             'fix_alpha': 0,
                             'alpha': 0.5,
                             'ncatG': 4})
        new_file = ctl.replace('.ctl', '.modify.ctl')
        with open(new_file, 'w') as f1:
            f1.write(new_text)
        params.append(
            f"cd {dirname(new_file)}; /home-user/thliao/software/paml4.9j/bin/codeml {basename(new_file)} > {basename(new_file).replace('.modify.ctl','.log')}")

    with mp.Pool(processes=30) as tp:
        _ = list(tqdm((tp.imap(run, params)), total=len(params)))

    # cat all rst2 into in.BV
    rst2_list = glob(join(tmp_indir, '*', 'rst2'))
    text = ''
    for rst2 in sorted(rst2_list,
                       key=lambda x: int(basename(x).replace('tmp', ''))):
        text += open(rst2).read()
    with open(join(odir, 'in.BV'), 'w') as f1:
        f1.write(text)


def final_mcmctree(template_ctl, inBV, odir, ndata):
    # for final mcmctree
    if not exists(odir):
        os.makedirs(odir)
    bd_paras = '1 1 0.1'
    rgene_gamma = '1 35 1'
    sigma2_gamma = '1 10 1'
    burnin = '2000'
    sampfreq = '2'
    nsample = '20000'
    seqfile_b = '.'+in_phyfile
    treefile_b = '.'+in_treefile
    outfile = './out'
    seqtype = 2
    clock = 2
    param = {'seqfile': seqfile_b,
             'treefile': treefile_b,
             'ndata': ndata,
             'seqtype': seqtype,
             'usedata': "2 in.BV 1",
             'outfile': outfile,
             'clock': clock,
             'BDparas': bd_paras,
             'rgene_gamma': rgene_gamma,
             'sigma2_gamma': sigma2_gamma,
             'burnin': burnin,
             'sampfreq': sampfreq,
             'nsample': nsample,
             'alpha': 0.5}
    text = modify(template_ctl,
                  **param)
    os.system(f'cp {inBV} {odir}')
    ofile = join(odir, '03_mcmctree.ctl')
    with open(ofile, 'w') as f1:
        f1.write(text)
    run(
        f"/home-user/thliao/software/paml4.9j/bin/mcmctree {ofile} > {ofile.replace('.ctl','.log')}")


def cli(tree, aln_dir, odir):
    pass


if __name__ == "__main__":
    cli()
