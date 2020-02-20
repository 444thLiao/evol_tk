import multiprocessing as mp
import os

import numpy as np
import pandas as pd
from tqdm import tqdm

from dating_workflow.step_script.dating_pro import run

target_ = ['set24', 'set13']
target_dir = './AR_set14'

for t in target_:
    # os.makedirs(f'AR_{t}',exist_ok=True)
    # os.makedirs(f'IR_{t}', exist_ok=True)
    # os.system(f'cp {target_dir}/04_mcmctree.ctl {target_dir}/in.BV ./AR_{t}')
    # os.system(f'cp {target_dir}/04_mcmctree.ctl {target_dir}/in.BV ./IR_{t}')
    cmd = f"""R -e "setwd('AR_{t}'); b = mcmc3r::make.beta(n=8, a=5, method='step-stones'); mcmc3r::make.bfctlf(b, ctlf='04_mcmctree.ctl', betaf='beta.txt')" """
    os.system(cmd)
    cmd = f"""R -e "setwd('IR_{t}'); b = mcmc3r::make.beta(n=8, a=5, method='step-stones'); mcmc3r::make.bfctlf(b, ctlf='04_mcmctree.ctl', betaf='beta.txt')" """
    os.system(cmd)

cmds = []
for t in target_:
    for _ in range(1, 9):
        cmds.append(f"cd AR_{t}/{_}/ ; mcmctree 04_mcmctree.ctl > run.log ")

with mp.Pool(processes=len(cmds)) as tp:
    r = list(tqdm(tp.imap(run, cmds), total=len(cmds)))


def get_v(rout):
    outs = rout.split('\n')
    idx1, idx2 = outs.index('$logml'), outs.index('$se')
    logml, se = map(lambda x: float(x.split(' ')[-1].strip()),
                    (outs[idx1 + 1], outs[idx2 + 1]))
    return logml, se


collect_df = pd.DataFrame(columns=['calibration set', 'model', 'Log marginal (s. d)', 'BF'])
target_ = ['set24', 'set13','set14',]
for t in target_:
    cmd = f"""R -e "setwd('AR_{t}'); AR<- mcmc3r::stepping.stones(); AR " """
    AR = os.popen(cmd).read()
    AR_logml, AR_se = get_v(AR)
    cmd = f"""R -e "setwd('IR_{t}'); IR<- mcmc3r::stepping.stones(); IR " """
    IR = os.popen(cmd).read()
    IR_logml, IR_se = get_v(IR)
    c = np.array([AR_logml, IR_logml])
    BF = np.exp(c - np.max(c))
    collect_df.loc[0, :] = [t, 'AR', f'{AR_logml} ({AR_se})', BF[0]]
    collect_df.loc[1, :] = [t, 'IR', f'{IR_logml} ({IR_se})', BF[1]]
