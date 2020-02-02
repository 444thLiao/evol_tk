import multiprocessing as mp
import os
import subprocess
from glob import glob
from os.path import *
from subprocess import check_call

import numpy as np
from tqdm import tqdm

# __file__ = '/home-user/thliao/script/evolution_relative/dating_workflow/step_script/get_subrate.py'
template_dir = abspath(join(dirname(dirname(__file__)), 'ctl_template'))
codeml_ctl = join(template_dir, 'codeml.ctl')
aaRatefile = join(template_dir, 'lg.dat')

paml_bin = "/home-user/software/paml/v4.9/paml4.9e/bin/codeml"


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
        # subprocess.check_output(cmd, shell=1)
        check_call(cmd,
                   shell=1,
                   stdout=open(log, 'w'))

    except subprocess.CalledProcessError as e:
        print('error', e.output)
    if log != '/dev/null':
        t = open(log, 'r', newline='\n').read().replace('\r', '\n')
        with open(log, 'w') as f1:
            f1.write(t)


odir = './'
aln_dir = '../cog25_single/187g_aln/'
genome_id = '../dating_for_187g.list'
in_treefile = abspath('./dating_for/phy_files/82g_subrate/187g_point.newick')

odir = './dating_for/phy_files/29cyano_subrate'
aln_dir = './cog25_single/29cyano_aln'
genome_id = '/home-user/thliao/data/cyano_basal/rawdata/assembly_ids.list'
in_treefile = abspath('./trees/final/29cyano.newick')

odir = './dating_for/phy_files/83g_subrate'
aln_dir = './cog25_single/83g_aln/'
genome_id = './dating_for_83g.list'
in_treefile = abspath('./dating_for/phy_files/83g_subrate/83g_point.newick')
if not exists(odir):
    os.makedirs(odir)

cmds = []
for f in glob(join(aln_dir, '*.trimal')):
    ofile = join(odir, basename(f).replace('.trimal', '.phy'))
    # if not exists(ofile):
    run(f'python3 ~/script/evolution_relative/dating_workflow/step_script/aln2phy.py -i {f} -o {ofile} -gl {genome_id}')
    if not exists(odir):
        os.makedirs(odir)

    seqfile_b = abspath(ofile)
    treefile_b = in_treefile
    outfile = f"./{basename(f).replace('.trimal', '.phy')}.out"
    seqtype = 2
    clock = 1
    aaRatefile = aaRatefile
    param = {'seqfile': seqfile_b,
             'treefile': treefile_b,
             'seqtype': seqtype,
             'outfile': outfile,
             'clock': clock,
             'aaRatefile': aaRatefile}
    text = modify(codeml_ctl, **param)
    ctl_f = join(odir, basename(f).replace('.trimal', '_codelml.ctl'))
    with open(ctl_f, 'w') as f1:
        f1.write(text)

    cmds.append(f"cd {dirname(ctl_f)}; {paml_bin} {basename(ctl_f)} ")

with mp.Pool(processes=30) as tp:
    r = list(tqdm(tp.imap(run, cmds), total=len(cmds)))

g2rate = {}
for outfile in glob(join(odir, '*.out')):
    gene_name = basename(outfile).partition('.')[0]
    rows = open(outfile).readlines()
    idx = [_ for _, r in enumerate(rows) if 'Substitution rate' in r]
    if not idx:
        continue
    idx = idx[0]
    rate = rows[idx + 1].strip()
    rate = float(rate.split('+-')[0].strip())
    g2rate[gene_name] = rate / 3.5 / 10
m = np.mean(list(g2rate.values()))
s = np.std(list(g2rate.values()))
import scipy.stats as stats

fit_alpha, fit_loc, fit_beta = stats.gamma.fit(list(g2rate.values()))
