import multiprocessing as mp
import os
from glob import glob
from os.path import *

# generate cal set
from tqdm import tqdm

from dating_workflow.step_script.dating_pro import modify, run

new_trees = glob('./dating_for/cal_tree/83g_set*.newick')
odir = './dating_for/83g/batch_prior'
cmds = []
for tree in new_trees:
    set_name = basename(tree).split('_')[-1].replace('.newick', '')
    prefix = basename(tree).split('_')[0]
    prefile = f'./dating_for/{prefix}/{prefix}_set1/prior/nodata_mcmctree.ctl'
    param = {
        # 'seqfile': seqfile_b,
        'treefile': abspath(tree),
        #          'ndata': ndata,
        #          'seqtype': seqtype,
        'usedata': "0",
        # 'outfile': './03_mcmctree.out',
        'clock': '2',
        #          'BDparas': bd_paras,
        'rgene_gamma': '1 100 1',
        #          'sigma2_gamma': sigma2_gamma,
        'burnin': 2000,
        'sampfreq': 20,
        'nsample': 20000,
        #          'alpha': 0.5
        # 'print':1
    }

    # modify these ctl
    onew_name = f'{prefix}_{set_name}_prior'
    if exists(prefile):
        text = modify(prefile, **param)
        final_odir = join(odir, onew_name)
        os.makedirs(final_odir, exist_ok=True)
        d = dirname(prefile)
        with open(join(final_odir, 'nodata_mcmctree.ctl'), 'w') as f1:
            f1.write(text)
        cmd = f"cd {join(final_odir)} ; mcmctree ./nodata_mcmctree.ctl > ./run.log "

        if not exists(join(final_odir, 'FigTree.tre')):
            cmds.append(cmd)

with mp.Pool(processes=len(cmds)) as tp:
    r = list(tqdm(tp.imap(run, cmds), total=len(cmds)))
