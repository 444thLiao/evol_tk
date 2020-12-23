"""
Get substitution rate
example. Still raw scripts
"""
import multiprocessing as mp
import os
import subprocess
from glob import glob
from os.path import *
from subprocess import check_call
import pandas as pd
import numpy as np
from tqdm import tqdm

# __file__ = '/home-user/thliao/script/evolution_relative/dating_workflow/step_script/get_subrate.py'
template_dir = abspath(join(dirname(dirname(__file__)), 'ctl_template'))
codeml_ctl = join(template_dir, 'codeml.ctl')
baseml_ctl = join(template_dir, 'baseml.ctl')
aaRatefile = join(template_dir, 'lg.dat')

paml_bin = "/home-user/software/paml/v4.9/paml4.9e/bin/codeml"
paml_bin_base = "/share/home-user/thliao/software/paml4.9e/bin/baseml"


phy_file = "/mnt/home-backup/thliao/plancto/update_calibrations/trees/phy_files/scheme1_cog25_nucl.phy"
tree_file = "/mnt/home-backup/thliao/plancto/update_calibrations/dating/cal_tree/78g_set10.newick"
param = {'seqfile': phy_file,
            'treefile': tree_file,
            'seqtype': 2,
            'outfile': './set10.out',
            'clock': 1,
            }
text = modify(codeml_ctl, **param)
with open('./set10_baseml.ctl', 'w') as f1:
        f1.write(text)
cmd = f"{paml_bin} ./set10_baseml.ctl"

cmds = []
for f in glob(join(aln_dir, '*.trimal')):
    seqfile_b = abspath(ofile)
    treefile_b = in_treefile
    outfile = f"./{basename(f).replace('.trimal', '.phy')}.out"
    seqtype = 2
    clock = 1
    aaRatefile = aaRatefile
    param = {'seqfile': phy_file,
             'treefile': tree_file,
             'seqtype': 2,
             'outfile': './set10.out',
             'clock': 1,
             #'aaRatefile': aaRatefile
             }
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


import plotly.graph_objs as go
import plotly.express as px
df1 = pd.DataFrame().from_dict({'rate':g2rate})

df = pd.concat([df1,df2],axis=0)
df.loc[:,'x'] = ['gene%s' % _ for _ in df.index.astype(str)]
fig = px.bar(df,x='x',y='rate',color='part',barmode='group')
fig.write_html("./compare_rate_nucl_prot.html", include_plotlyjs='cdn')