from glob import glob
import multiprocessing as mp
import os
from os.path import *
from dating_workflow.step_script.dating_pro import modify,run,final_mcmctree
from tqdm import tqdm

in_pattern = './dating_for/198g*'



process_dirs = []
for indir in glob(in_pattern):
    if isdir(indir):
        process_dirs.append(indir)

param = {
    # 'seqfile': seqfile_b,
    #          'treefile': treefile_b,
    #          'ndata': ndata,
    #          'seqtype': seqtype,
    #          'usedata': "2 in.BV 1",
    #          'outfile': outfile,
    #          'clock': clock,
    #          'BDparas': bd_paras,
    #          'rgene_gamma': rgene_gamma,
    #          'sigma2_gamma': sigma2_gamma,
    #          'burnin': burnin,
             'sampfreq': 5,
    #          'nsample': nsample,
    #          'alpha': 0.5
             }

# modify these ctl
for d in process_dirs:
    ori_ctl = join(d,'mcmc_for','03_mcmctree.ctl')
    if exists(ori_ctl):
        text = modify(ori_ctl,**param)
        with open(join(d,'mcmc_for','04_mcmctree.ctl'),'w') as f1:
            f1.write(text)

cmds = []
for d in process_dirs:
    cmd = f"cd {join(d,'mcmc_for')} ; mcmctree ./04_mcmctree.ctl"
    cmds.append(cmd)
    
    
with mp.Pool(processes=10) as tp:
    r = list(tqdm(tp.imap(run, cmds), total=len(cmds)))