import multiprocessing as mp
import os
from os.path import exists,join
from subprocess import check_call
import click
import numpy as np
import pandas as pd
from tqdm import tqdm

from dating_workflow.bin.dating_pro import run, modify
from glob import glob

def run_mcmc3r(indir,odir,intree,clock=2,usedata=1):
    if not exists(odir):
        os.system(f'mkdir -p {odir}')
    ctl = glob(f"{indir}/*.ctl")
    if not ctl:
        raise IOError('not ctl file')
    
    os.system(f'cp {ctl[0]} {odir}/04_mcmctree.ctl')
    os.system(f'ln -sf `realpath {indir}/in.BV` {odir}/')
    param = {'treefile': intree,
              'clock': clock,
              'usedata': '1' if usedata ==1 else '2 ../in.BV 1',
            #  'seqtype': '0'
                 }
    text = modify(f'{odir}/04_mcmctree.ctl',
                    **param)
    with open(f'{odir}/04_mcmctree.ctl', 'w') as f1:
        f1.write(text)
    cmd = f"""/home-user/thliao/anaconda3/envs/r_env/bin/R -e "setwd('{odir}'); b = mcmc3r::make.beta(n=8, a=5, method='step-stones'); mcmc3r::make.bfctlf(b, ctlf='04_mcmctree.ctl', betaf='beta.txt')" """
    check_call(cmd,shell=1,executable='/home-user/thliao/anaconda3/bin/zsh')

    _ctl = "04_mcmctree.ctl"
    cmds = []
    for _ in range(1, 9):
        cmds.append(f"cd {odir}/{_}/ ; mcmctree {_ctl} > run.log ")
    return cmds

def main():
    pass


if __name__ == '__main__':
    # target_ = ['set24', 'set13', 'set1', 'set14', 'set25']
    target_ = ['set33', 'set34', 'set35', 'set36', 'set37']
    target_dir = './AR_set1'

    for t in target_:
        for model in ['IR', 'AR']:
            if not exists(f"./{model}_{t}"):
                os.makedirs(f"./{model}_{t}")
            os.system(f'cp {target_dir}/03_mcmctree.ctl ./{model}_{t}/')
            param = {'treefile': "/share/home-user/thliao/data/plancto/dating_for/cal_tree/83g_set1.newick".replace('set1', t),
                     'clock': 2 if model == 'IR' else 3}
            text = modify(f'./{model}_{t}/03_mcmctree.ctl',
                          **param)
            with open(f'./{model}_{t}/03_mcmctree.ctl', 'w') as f1:
                f1.write(text)
            cmd = f"""R -e "setwd('{model}_{t}'); b = mcmc3r::make.beta(n=8, a=5, method='step-stones'); mcmc3r::make.bfctlf(b, ctlf='03_mcmctree.ctl', betaf='beta.txt')" """
            os.system(cmd)
        cmd = f"""R -e "setwd('IR_{t}'); b = mcmc3r::make.beta(n=8, a=5, method='step-stones'); mcmc3r::make.bfctlf(b, ctlf='03_mcmctree.ctl', betaf='beta.txt')" """
        os.system(cmd)

    _ctl = "03_mcmctree.ctl"
    cmds = []
    for t in target_:
        for _ in range(1, 9):
            for model in ['IR', 'AR']:
                cmds.append(f"cd {model}_{t}/{_}/ ; mcmctree {_ctl} > run.log ")

    with mp.Pool(processes=len(cmds)) as tp:
        r = list(tqdm(tp.imap(run, cmds), total=len(cmds)))


    def get_v(rout):
        outs = rout.split('\n')
        idx1, idx2 = outs.index('$logml'), outs.index('$se')
        logml, se = map(lambda x: float(x.split(' ')[-1].strip()),
                        (outs[idx1 + 1], outs[idx2 + 1]))
        return logml, se


    collect_df = pd.DataFrame(columns=['calibration set', 'model', 'Log marginal (s. d)', 'BF'])
    count = 0
    # for t in ['AR_sim','IR_sim']:
    cmd = f"""/home-user/thliao/anaconda3/envs/r_env/bin/R -e "setwd('./dating/sys_testing/model_selection/C9_no_euk/AR_sim'); AR<- mcmc3r::stepping.stones(); AR " """
    AR = os.popen(cmd).read()
    AR_logml, AR_se = get_v(AR)
    cmd = f"""/home-user/thliao/anaconda3/envs/r_env/bin/R -e "setwd('./dating/sys_testing/model_selection/C9_no_euk/IR_sim'); IR<- mcmc3r::stepping.stones(); IR " """
    IR = os.popen(cmd).read()
    IR_logml, IR_se = get_v(IR)
    c = np.array([AR_logml, IR_logml])
    BF = np.exp(c - np.max(c))
    collect_df.loc[count, :] = [t, 'AR', f'{AR_logml} ({AR_se})', BF[0]]
    collect_df.loc[count + 1, :] = [t, 'IR', f'{IR_logml} ({IR_se})', BF[1]]
    count += 2
