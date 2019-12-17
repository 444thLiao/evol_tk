import multiprocessing as mp
import subprocess
from glob import glob
from os.path import *
from subprocess import check_call

from tqdm import tqdm

# TODO: run each locus indepently
# test is it the fault of which locus?
indir = expanduser("~/data/nitrification_for/dating_for/cog25_single/243g")

all_locus = glob(join(indir, '*.aln'))
all_locus = [basename(_).replace('.aln', '') for _ in all_locus]

# print(all_locus)
for each in all_locus:
    ofile = f"./locus_each/normal_rmI_25g_nofill_{each}.trimal"
    cmd = f"python3 ~/script/evolution_relative/dating_workflow/toolkit/concat_aln.py -i ~/data/nitrification_for/dating_for/cog25_single/243g -ct phy -gl ~/data/nitrification_for/dating_for/bac120_annoate/remained_ids_fullv1.list -o {ofile} -s trimal -no_graph -no_fill -rm_I -genel {each}"
    if not exists(ofile):
        check_call(cmd, shell=1)


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


params = []
for each in all_locus:
    ofile = f"./locus_each/normal_rmI_25g_nofill_{each}.phy"
    cmd = f"python3 ~/script/evolution_relative/dating_workflow/step_script/dating_pro.py -i {ofile} -it ./243g_120gene_5cal.newick -o ./locus_each/normal_rmI_25g_nofill_{each} -only_prior"
    params.append(cmd)

with mp.Pool(processes=30) as tp:
    _ = list(tqdm((tp.imap(run, params)), total=len(params)))
