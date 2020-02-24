import multiprocessing as mp
import os
from os.path import exists

import numpy as np
import pandas as pd
from tqdm import tqdm

from dating_workflow.step_script.dating_pro import run, modify

target_ = ['set24', 'set13', 'set1', 'set14', 'set25']

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
for t in target_:
    cmd = f"""R -e "setwd('AR_{t}'); AR<- mcmc3r::stepping.stones(); AR " """
    AR = os.popen(cmd).read()
    AR_logml, AR_se = get_v(AR)
    cmd = f"""R -e "setwd('IR_{t}'); IR<- mcmc3r::stepping.stones(); IR " """
    IR = os.popen(cmd).read()
    IR_logml, IR_se = get_v(IR)
    c = np.array([AR_logml, IR_logml])
    BF = np.exp(c - np.max(c))
    collect_df.loc[count, :] = [t, 'AR', f'{AR_logml} ({AR_se})', BF[0]]
    collect_df.loc[count + 1, :] = [t, 'IR', f'{IR_logml} ({IR_se})', BF[1]]
    count += 2
