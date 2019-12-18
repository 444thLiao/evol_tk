"""
get result from multistate+MCMC
summarizing into iTOL annotation files.
"""
import io
from collections import defaultdict

import pandas as pd

from api_tools.itol_func import pie_chart


# infile = 'm2nm.txt.Log.txt'

def summaized_r(complex_f, simple_f, key=''):
    complex_df = get_df(complex_f, key=key)
    simple_df = get_df(simple_f, key=key)

    complex_lh = complex_df.iloc[-1, 1]
    simple_lh = simple_df.iloc[-1, 1]

    lrt = 2 * (complex_lh - simple_lh)
    text = f"Log BF = 2({complex_lh}- {simple_lh})\nLog BF = {lrt}"
    return text


def get_df(infile, key='Iteration'):
    rows = open(infile).readlines()
    header = [(idx, _) for idx, _ in enumerate(rows) if _.startswith(key)][-1]
    start_at = header[0]
    headers = header[1].strip('\n')
    headers = headers.split('\t')

    result_df = pd.read_csv(io.StringIO(''.join(rows[start_at:])), sep='\t')
    return result_df


def get_result(infile,cat2info={}):
    result_df = get_df(infile)

    mean_vals = result_df.mean()
    mean_v = mean_vals.to_dict()

    transition_rate = {k: v
                       for k, v in mean_v.items() if k.startswith('q')}

    n2cat2prob = defaultdict(lambda: defaultdict(dict))
    for key, v in mean_v.items():
        if ' P(' in key:
            node_name = key.split(' ')[0]
            if node_name == 'Root':
                node_name = 'OROOT'  # for itol
            cat = key.split('(')[-1].strip(')')
            n2cat2prob[node_name][cat] = v

    cat2info = {"M": '#0000ff',
                "N": '#D68529'}

    text = pie_chart(n2cat2prob, cat2info, )

    return text
