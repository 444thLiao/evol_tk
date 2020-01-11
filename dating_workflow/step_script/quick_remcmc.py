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
              'usedata': "2 ../in.BV 1",
              'outfile': './04_mcmctree.out',
    #          'clock': clock,
    #          'BDparas': bd_paras,
              'rgene_gamma': '1 100 1',
    #          'sigma2_gamma': sigma2_gamma,
    #          'burnin': burnin,
             'sampfreq': 10,
    #          'nsample': nsample,
    #          'alpha': 0.5
             }

# modify these ctl
onew_name = 'repeat_rg1_100'
for d in process_dirs:
    ori_ctl = join(d,'mcmc_for','03_mcmctree.ctl')
    if exists(ori_ctl):
        text = modify(ori_ctl,**param)
        os.makedirs(join(d,'mcmc_for',onew_name),exist_ok=True)
        with open(join(d,'mcmc_for',onew_name,'04_mcmctree.ctl'),'w') as f1:
            f1.write(text)
# remove previous
# for d in process_dirs:
#     ori_f = join(d,'mcmc_for','mcmc.txt')
#     if exists(ori_f):
#         os.system(f'rm -rf {ori_f}',shell=1)
#     ori_f = join(d,'mcmc_for','FigTree.tre')
#     if exists(ori_f):
#         os.system(f'rm -rf {ori_f}',shell=1)
#     ori_f = join(d,'mcmc_for','SeedUsed')
#     if exists(ori_f):
#         os.system(f'rm -rf {ori_f}',shell=1)        
        
cmds = []
for d in process_dirs:
    cmd = f"cd {join(d,'mcmc_for',onew_name,)} ; mcmctree ./04_mcmctree.ctl"
    cmds.append(cmd)
    
    
with mp.Pool(processes=10) as tp:
    r = list(tqdm(tp.imap(run, cmds), total=len(cmds)))