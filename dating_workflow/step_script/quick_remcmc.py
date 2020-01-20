import multiprocessing as mp
import os
from glob import glob
from os.path import *

# generate cal set
import pandas as pd
from tqdm import tqdm

from api_tools.for_tree.format_tree import add_cal_api
from dating_workflow.step_script.dating_pro import modify, run

a = pd.read_excel('calibrations_sets.xlsx', index_col=0)
mapping_dict = {
    'Node 1': 'GCA_000011385.1|GCA_002746535.1',
    'Node 2': 'GCA_000011385.1|GCA_000013205.1',
    'Node 3': 'GCA_000196515.1|GCA_001548455.1',
    'Node 4': 'GCA_000317575.1|GCA_000317025.1', }
name_dict = {
    'Node 1': 'root',
    'Node 2': 'gleobacteria',
    'Node 3': 'Nostocales',
    'Node 4': 'pleurocapsales',
}

for s, row in a.iterrows():
    rows = []
    for node, t in row.items():
        if not pd.isna(t):
            rows.append(f"{mapping_dict[node]}\t{str(t)}\t{name_dict[node]}")
    with open(f'./cal_{s}.txt', 'w') as f1:
        f1.write('\n'.join(rows))

redo = False
# batch cal set
ori_newick = './trees/final/198g_merged.newick'
new_trees = []
for cal_set in glob('./dating_for/calibrations_set/cal_set*.txt'):
    set_name = basename(cal_set).split('_')[-1].replace('.txt', '')
    if redo:
        add_cal_api(ori_newick,
                    f'./dating_for/cal_tree/198g_{set_name}.newick',
                    cal_set,
                    format=3)
    new_trees.append(abspath(f'./dating_for/cal_tree/198g_{set_name}.newick'))

ori_newick = './trees/final/187g_merged.newick'
for cal_set in glob('./dating_for/calibrations_set/cal_set*.txt'):
    set_name = basename(cal_set).split('_')[-1].replace('.txt', '')
    if redo:
        add_cal_api(ori_newick,
                f'./dating_for/cal_tree/187g_{set_name}.newick',
                cal_set,
                format=3)
    new_trees.append(abspath(f'./dating_for/cal_tree/187g_{set_name}.newick'))

ori_newick = './trees/final/184g_merged.newick'
for cal_set in glob('./dating_for/calibrations_set/cal_set*.txt'):
    set_name = basename(cal_set).split('_')[-1].replace('.txt', '')
    if redo:
        add_cal_api(ori_newick,
                f'./dating_for/cal_tree/184g_{set_name}.newick',
                cal_set,
                format=3)
    new_trees.append(abspath(f'./dating_for/cal_tree/184g_{set_name}.newick'))

odir = './dating_for/modify_rg'
cmds = []
for tree in new_trees:
    set_name = basename(tree).split('_')[-1].replace('.newick', '')
    prefix = basename(tree).split('_')[0]
    prefile = f'./dating_for/{prefix}/{prefix}_3cal_set1/mcmc_for/03_mcmctree.ctl'
    param = {
        # 'seqfile': seqfile_b,
        'treefile': abspath(tree),
        #          'ndata': ndata,
        #          'seqtype': seqtype,
        'usedata': "2 in.BV 1",
        'outfile': './04_mcmctree.out',
        #          'clock': clock,
        #          'BDparas': bd_paras,
        'rgene_gamma': '13 1246.5 1',
        #          'sigma2_gamma': sigma2_gamma,
        'burnin': 2000,
        'sampfreq': 10,
        'nsample': 20000,
        #          'alpha': 0.5
    }

    # modify these ctl
    for repeat_n in ['run1', 'run2']:
        onew_name = f'repeat_{prefix}_{set_name}_{repeat_n}'
        if exists(prefile):
            text = modify(prefile, **param)
            final_odir = join(odir, onew_name)
            os.makedirs(final_odir, exist_ok=True)
            d = dirname(prefile)
            os.system(f"ln -s {abspath(join(d, 'in.BV'))} {final_odir}/ ")
            with open(join(final_odir, '04_mcmctree.ctl'), 'w') as f1:
                f1.write(text)
            cmd = f"cd {join(final_odir)} ; mcmctree ./04_mcmctree.ctl"
            cmds.append(cmd)

# remove previous
for d in cmds:
    d = d.split(' ')[1]
    ori_f = join(d,'mcmc.txt')
    if exists(ori_f):
        os.system(f'rm -rf {ori_f}')
    ori_f = join(d,'FigTree.tre')
    if exists(ori_f):
        os.system(f'rm -rf {ori_f}')
    ori_f = join(d,'SeedUsed')
    if exists(ori_f):
        os.system(f'rm -rf {ori_f}')
    ori_f = join(d, '*.out')
    os.system(f'rm -rf {ori_f}')


with mp.Pool(processes=len(cmds)) as tp:
    r = list(tqdm(tp.imap(run, cmds), total=len(cmds)))
