"""
read the generated output of mcmctree

"""
import pandas as pd
def read_outfile(f):
    f = open(f).read().split('\n')
    head = 'Posterior mean (95% Equal-tail CI) (95% HPD CI) HPD-CI-width'
    if head not in f:
        print('no header, maybe not finish')
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
        if row:
            vals = row.split(' ')
            vals = [_ for _ in vals if _ and _ not in '(),']
            posterior_mean, equal_tailed_95p_down, equal_tailed_95p_up = map(format_v, vals[1:4])
            CI_collect.append((equal_tailed_95p_up - equal_tailed_95p_down))
            mean_collect.append(posterior_mean)
            idx.append(vals[0])
            CIs.append('%s - %s' % (equal_tailed_95p_down,equal_tailed_95p_up))
    df = pd.DataFrame()
    df.loc[:, 'CI_width'] = CI_collect
    df.loc[:,'CIs'] = CIs
    df.loc[:, 'Posterior mean time (100 Ma)'] = mean_collect
    df.index = idx
    return df