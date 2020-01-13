from glob import glob
import multiprocessing as mp
import os
from os.path import *
from dating_workflow.step_script.dating_pro import modify,run,final_mcmctree
from tqdm import tqdm
# batch cal set
ori_newick = './trees/final/198g_merged.newick'
new_trees = []
for cal_set in glob('./dating_for/calibrations_set/cal_set*.txt'):
    set_name = basename(cal_set).split('_')[-1].replace('.txt','')
    add_cal_api(ori_newick, 
                f'./dating_for/cal_tree/198g_{set_name}.newick', 
                cal_set, 
                format=3)
    new_trees.append(abspath(f'./dating_for/cal_tree/198g_{set_name}.newick'))
    
ori_newick = './trees/final/187g_merged.newick'
for cal_set in glob('./dating_for/calibrations_set/cal_set*.txt'):
    set_name = basename(cal_set).split('_')[-1].replace('.txt','')
    add_cal_api(ori_newick, 
                f'./dating_for/cal_tree/187g_{set_name}.newick', 
                cal_set, 
                format=3)
    new_trees.append(abspath(f'./dating_for/cal_tree/187g_{set_name}.newick'))

in_pattern = './dating_for/187g'
process_dirs = []
for indir in glob(in_pattern):
    if isdir(indir):
        process_dirs.append(indir)
        
process_dirs = ['./dating_for/187g']
cmds = []
for tree in new_trees:
    set_name = basename(tree).split('_')[-1].replace('.newick','')
    param = {
        # 'seqfile': seqfile_b,
                 'treefile': tree,
        #          'ndata': ndata,
        #          'seqtype': seqtype,
                'usedata': "2 ../in.BV 1",
                'outfile': './04_mcmctree.out',
        #          'clock': clock,
        #          'BDparas': bd_paras,
                'rgene_gamma': '1 100 1',
        #          'sigma2_gamma': sigma2_gamma,
                'burnin': 50000,
                'sampfreq': 25,
        #          'nsample': nsample,
        #          'alpha': 0.5
                }

    # modify these ctl
    for repeat_n in ['run1','run2']:
        onew_name = f'repeat_{set_name}_{repeat_n}'
        for d in process_dirs:
            ori_ctl = join(d,'mcmc_for','03_mcmctree.ctl')
            if exists(ori_ctl):
                text = modify(ori_ctl,**param)
                os.makedirs(join(d,'mcmc_for',onew_name),exist_ok=True)
                with open(join(d,'mcmc_for',onew_name,'04_mcmctree.ctl'),'w') as f1:
                    f1.write(text)
                cmd = f"cd {join(d,'mcmc_for',onew_name,)} ; mcmctree ./04_mcmctree.ctl"
                cmds.append(cmd)
                
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
        
    
with mp.Pool(processes=20) as tp:
    r = list(tqdm(tp.imap(run, cmds), total=len(cmds)))