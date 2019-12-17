import os
from collections import Counter
from glob import glob
from os.path import join, basename, exists, dirname

import click
import pandas as pd


def get_shared_ko(pattern, outfile, num, substite_org_name, indir):
    if '/' not in outfile:
        outfile = './' + outfile
    file_list = glob(str(join(indir, pattern)))
    if not file_list:
        file_list = [f for f in file_list if basename(f) in pattern.split(',')]
    if not exists(dirname(outfile)):
        os.makedirs(dirname(outfile))
    if not file_list:
        exit(f"No file detected according your input parameters 'indir/{pattern}', please check it")
    collect_dfs = []
    for f in file_list:
        try:
            _df = pd.read_csv(f, sep='\t', header=None, index_col=0)
        except pd.errors.EmptyDataError:
            print(f"empty file {f}")
            continue
        _df.loc[:, 'name'] = basename(f)
        collect_dfs.append(_df)

    counter_set = Counter(
        [ko
         for each_df in collect_dfs
         for ko in each_df.index
         # if each_df.loc[ko,1]!=0 or each_df.loc[ko,2]!=0
         if each_df.loc[ko, 3] != 0
         # uncomment above line could filter out the zero gain/loss rows
         ])
    counter_set = dict(counter_set)
    add_set = {ko: 0 for each_df in collect_dfs for ko in each_df.index if ko not in counter_set}
    counter_set.update(add_set)
    if num is None:
        nums = list(set(counter_set.values()))
    else:
        nums = [num]

    collect_df_list = []
    for num in nums:
        num = int(num)
        shared_ko_list = [k for k, v in counter_set.items()
                          if v == num]
        new_collect_dfs = [sub_df.reindex(shared_ko_list) for sub_df in collect_dfs]
        final_df = pd.concat(new_collect_dfs, axis=0)
        final_df = final_df.loc[~final_df.isna().all(1), :]
        if final_df.shape[0] == 0:
            continue
        # final_df.loc[:,'appear times'] = ''
        final_df.loc[:, 'appear times'] = num
        collect_df_list.append(final_df)
    collect_df_list = [_ for _ in collect_df_list if _.shape[0] != 0]
    if len(collect_df_list) == 0:
        print('No valid dataframe found ,exit......')
        return
    final_df = pd.concat(collect_df_list, axis=0)
    final_df = final_df.sort_values('appear times', ascending=False)
    if substite_org_name:
        rows = open(substite_org_name, 'r').read().split('\n')
        rows = [_ for _ in rows if _]
        org2new_name = dict([row.split('\t') for row in rows])
        final_df.loc[:, 'name'] = [org2new_name.get(_, 'Unknown') for _ in final_df.loc[:, 'name']]
    if '.xlsx' in outfile or '.xls' in outfile:
        final_df.to_excel(outfile, index=1)
    elif '.tab' in outfile or '.tsv' in outfile:
        final_df.to_csv(outfile, index=1, sep=',')
    else:
        final_df.to_csv(outfile, index=1, sep=',')


#    return shared_ko_list


@click.command()
@click.option('-p', 'pattern', help='regular expression or comma separated file name')
@click.option('-o', 'outfile', help='path of output file')
@click.option('-n', 'number', help='number of times you want to count', default=None)
@click.option('-i', 'indir', required=False, default='./', help='input directory of your pattern')
@click.option('-sub', 'substite_org_name', help='accept a file to replace the org name', default=None)
def cli(pattern, outfile, number, indir, substite_org_name):
    get_shared_ko(pattern, outfile, number, substite_org_name, indir=indir, )


if __name__ == "__main__":
    cli()
