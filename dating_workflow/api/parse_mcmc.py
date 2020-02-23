from glob import glob
from os.path import *

import click
import numpy as np
import pandas as pd
from ete3 import Tree
from tqdm import tqdm

from dating_workflow.debug_for.draw_infinite_plot import get_plot
import warnings
warnings.filterwarnings('ignore')

def cal_ESS(df):
    N = int(df.shape[0] * 0.1)
    df = df.iloc[N:, :]
    N = df.shape[0]

    def f(t1):
        def cal(x, k):
            return np.corrcoef(x[k:], x[:-k])[0][1]

        V = np.array([cal(t1.values, k) for k in range(1, N - 1)])
        ess = N / (1 + sum(V))
        return ess

    ess_array = df.apply(f, axis=0)


def get_CI(f):
    """

    :param f: log file?
    :return:
    """
    f = open(f).read().split('\n')
    head = 'Posterior means (95% Equal-tail CI) (95% HPD CI) HPD-CI-width'
    if head not in f:
        return None
    idx = f.index(head)
    remained_txt = f[idx + 1:]

    def format_v(x):
        x = x.strip('(')
        x = x.strip(')')
        x = x.strip(',')
        return float(x)

    idx = []
    CIs = []
    mean_collect = []
    CI_collect = []
    for row in remained_txt:
        if row.startswith('t_n') or row.startswith('lnL'):
            vals = row.split(' ')
            vals = [_ for _ in vals if _ and _ not in '(),']
            posterior_mean, equal_tailed_95p_down, equal_tailed_95p_up = map(format_v, vals[1:4])
            CI_collect.append((equal_tailed_95p_up - equal_tailed_95p_down))
            mean_collect.append(posterior_mean)
            idx.append(vals[0])
            CIs.append('%s - %s' % (equal_tailed_95p_down, equal_tailed_95p_up))
    df = pd.DataFrame()
    df.loc[:, 'CI_width'] = CI_collect
    df.loc[:, 'CIs'] = CIs
    df.loc[:, 'Posterior mean time (100 Ma)'] = mean_collect
    df.index = idx
    return df


def get_node_name(f):
    matched_row = ''
    match = False
    for row in open(f):
        row = row.strip('\n').strip()
        if match and not row:
            break
        if row.startswith('Species tree for FigTree.  Branch lengths = posterior mean times; 95% CIs = labels'):
            match = True
            continue
        if match:
            matched_row = row
    t = Tree(matched_row.replace(' ', ''), format=8)
    for l in t.get_leaves():
        l.name = l.name.partition('_')[-1]
    return t


def main(indir, ns, groupname, odir):
    tmp_df = pd.DataFrame()
    collect_ = {}

    processed_dir = list(glob(join(indir, '*run1')))
    name = ''
    for each_dir in tqdm(processed_dir):
        outfile = glob(join(each_dir, '*.out'))
        if not outfile:
            continue
        else:
            outfile = outfile[0]
        mcmc = join(each_dir, 'mcmc.txt')
        log = join(each_dir, 'run.log')
        if not exists(log):
            log = glob(join(each_dir, '*.log'))[0]
        if exists(join(each_dir, 'FigTree.tre')):
            if not name:
                t = get_node_name(outfile)
                name = 't_n%s' % t.get_common_ancestor(ns).name
            df = pd.read_csv(mcmc, sep='\t', index_col=0)
            set_name = basename(dirname(outfile)).partition('_')[-1]
            if not set_name:
                set_name = basename(dirname(outfile))
            collect_[set_name] = df.loc[:, name]
            df = get_CI(log)
            tmp_df.loc[set_name, groupname] = '%s (%s) ' % (df.loc[name, 'Posterior mean time (100 Ma)'],
                                                            df.loc[name, 'CIs'])
            tmp_df.loc[set_name, 'ROOT'] = '%s (%s) ' % (df.loc[:, 'Posterior mean time (100 Ma)'].values[0],
                                                         df.loc[:, 'CIs'].values[0])
            tmp_df.loc[set_name, 'lnL'] = '%s (%s) ' % (df.loc["lnL", 'Posterior mean time (100 Ma)'],
                                                        df.loc["lnL", 'CIs'])

    tmp_df.loc[:, 'num_set'] = [int(_.split('_')[1].replace('set', ''))
                                if 'run' in _ else 0
                                for _ in tmp_df.index]
    tmp_df = tmp_df.sort_values('num_set')
    # tmp_df.columns = ["divergence time/100Mya (CI)"]

    odir = join(odir, 'parsed_mcmc_result')  # './dating_for/83g/clock2_infinite_plot'
    pattern = join(indir, '*_run1', 'run.log')  # "./dating_for/83g/clock2_diff_cal/*_run1/run.log"
    df = get_plot(pattern, odir)

    writer = pd.ExcelWriter(join(odir, 'mcmc.xlsx'), engine='xlsxwriter')
    tmp_df.to_excel(writer, index_label='Calibration sets', sheet_name='Posterior time')
    df.to_excel(writer, index_label='Calibration sets', sheet_name='infinite site plots')
    writer.save()


@click.command()
@click.option("-i", 'indir', help='dir have multiple calibration set')
@click.option("-ns", 'targe_group', default=None)
@click.option('-o', 'odir')
def cli(indir):
    pass


if __name__ == '__main__':
    cli()

    ns = ['GCA_001828545.1', 'GCA_004282745.1']
    indir = './dating_for/83g/clock2_diff_cal/'
    main(indir=indir, ns=ns, groupname='Anammox group', odir=indir)
