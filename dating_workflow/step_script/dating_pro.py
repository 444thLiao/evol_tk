"""
20191202 Liao Tianhua
This script is mainly wrriten for debugging the mcmctree problems.
Of course it could run mainly mcmctree from input tree, alignment, aaRatemodel.

For now, it is designed for protein alignment only. 

For better understanding the usage of this script, see usage_dating.md
"""

from ete3 import Tree
import click
from glob import glob
import multiprocessing as mp
from subprocess import check_call, check_output
import subprocess
import os
from os.path import *
from tqdm import tqdm
import click

# __file__ = '/home-user/thliao/script/evolution_relative/dating_workflow/step_script/dating_pro.py'
# template file dir
template_dir = join(dirname(dirname(__file__)),'ctl_template')
mcmc_ctl = join(template_dir,'mcmctree.ctl')
codeml_ctl = join(template_dir,'codeml.ctl')
aaRatefile = join(template_dir,'lg.dat')

paml_bin = "/home-user/thliao/software/paml4.9j/bin"
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
        log = '/dev/null'
    else:
        cmd, log = args
    try:
        #subprocess.check_output(cmd, shell=1)
        check_call(cmd,
                   shell=1,
                   stdout=open(log,'w'))

    except subprocess.CalledProcessError as e:
        print('error',e.output)
    if log != '/dev/null':
        t = open(log,'r',newline='\n').read().replace('\r','\n')
        with open(log,'w') as f1:
            f1.write(t)

def get_num_phy_file(in_phyfile):
    ndata = 0
    for row in open(in_phyfile):
        row = row.strip('\n').strip().split(' ')
        row = [_ for _ in row if _]
        if all([_.isnumeric() for _ in row]):
            ndata +=1
    return ndata

def generate_tmp(in_phyfile,in_treefile,odir,ndata,template_ctl=mcmc_ctl):
    # template_ctl_01 = './01_mcmctree.ctl'
    if not exists(odir):
        os.makedirs(odir)
    new_01_ctl = join(odir, '01_mcmctree_modify.ctl')
    params = {'ndata': ndata,
              'seqfile': in_phyfile,
              'treefile': in_treefile,
              'outfile': './01_out'}
    text = modify(template_ctl, **params)
    with open(new_01_ctl, 'w') as f1:
        f1.write(text)
    run(f"export PATH=''; cd {odir}; {paml_bin}/mcmctree 01_mcmctree_modify.ctl")


def collecting_tmp(tmp_indir, odir):
    if not exists(odir):
        os.makedirs(odir)
    ctls = glob(join(tmp_indir, 'tmp0*.ctl'))
    for ctl in ctls:
        name = basename(ctl).replace('.ctl', '')
        os.makedirs(join(odir, f'{name}'), exist_ok=1)
        os.system(f'mv {tmp_indir}/{name}.* {odir}/{name}/')


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
        params.append((
            f"cd {dirname(new_file)}; {paml_bin}/codeml {basename(new_file)}", new_file.replace('.modify.ctl','.log')))

    with mp.Pool(processes=30) as tp:
        _ = list(tqdm((tp.imap(run, params)), total=len(params)))

    # cat all rst2 into in.BV
    rst2_list = glob(join(tmp_indir, '*', 'rst2'))
    text = ''
    for rst2 in sorted(rst2_list,
                       key=lambda x: int(basename(dirname(x)).replace('tmp', ''))):
        text += open(rst2).read()
    with open(join(odir, 'in.BV'), 'w') as f1:
        f1.write(text)


def final_mcmctree(inBV,in_phyfile,in_treefile, odir, ndata,template_ctl=mcmc_ctl):
    # for final mcmctree
    if not exists(odir):
        os.makedirs(odir)
    bd_paras = '1 1 0.1'
    rgene_gamma = '1 35 1'
    sigma2_gamma = '1 10 1'
    burnin = '2000'
    sampfreq = '2'
    nsample = '20000'
    seqfile_b = in_phyfile
    treefile_b = in_treefile
    outfile = './03_mcmctree.out'
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
    if not exists(f'{odir}/in.BV'):
        os.system(f'cp {inBV} {odir}')
    ofile = join(odir, '03_mcmctree.ctl')
    with open(ofile, 'w') as f1:
        f1.write(text)
    run( (f"cd {dirname(ofile)}; {paml_bin}/mcmctree 03_mcmctree.ctl ",
          ofile.replace('.ctl','.log'))  )


def main(in_phyfile, in_treefile, total_odir,run_tmp=True):
    if not exists(total_odir):
        os.makedirs(total_odir)
    ndata = get_num_phy_file(in_phyfile)
    mcmc_for_dir = join(total_odir,'mcmc_for')
    tmp_odir = join(total_odir,'tmp_files')
    if run_tmp:
        
    
        generate_tmp(in_phyfile, 
                    in_treefile,
                    tmp_odir,
                    ndata)
        collecting_tmp(tmp_odir,
                    tmp_odir)
        
        run_each_tmp(tmp_odir,
                    mcmc_for_dir)
    final_mcmctree(inBV=join(mcmc_for_dir,'in.BV'),
                    in_phyfile=in_phyfile,
                    in_treefile=in_treefile,
                    odir=mcmc_for_dir,
                    ndata=ndata)

def process_path(path):
    if not '/' in path:
        path = './' + path
    path = abspath(expanduser(path))
    return path


@click.command()
@click.option('-i','--in_phy','in_phyfile')
@click.option('-it','--in_tree','in_treefile')
@click.option('-o','odir')
@click.option('-no_tmp','run_tmp',is_flag=True, default=True)
def cli(in_phyfile, in_treefile, odir,run_tmp):
    in_phyfile = process_path(in_phyfile)
    in_treefile = process_path(in_treefile)
    main(in_phyfile, in_treefile, total_odir=odir,run_tmp=run_tmp)


if __name__ == "__main__":
    cli()
