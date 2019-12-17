import sys
from collections import defaultdict

import pandas as pd
from tqdm import tqdm

assert len(sys.argv) == 3

input_xlsx = sys.argv[1]
output_xlsx = sys.argv[2]

ori_df = pd.read_excel(input_xlsx, sheet_name=0)

['appear\ntimes', 'database', 'KO number', 'extant', 'ancestral',
 'gain or loss', 'genename', 'gene function description', 'nod to nif',
 'pathway', 'pathway description']

names = ori_df.index[(ori_df.loc[:, 'gain or loss'].isna() &
                      ~ori_df.loc[:, 'appear\ntimes'].isna())]
group_idx = defaultdict(list)
group_name = ''
for _, row in ori_df.iterrows():
    if _ not in names and not row.isna().all():
        group_idx[group_name].append(_)
    elif _ in names:
        group_name = ori_df.loc[_, 'appear\ntimes']
    else:
        pass

if not group_idx:
    pass

result_df = pd.DataFrame(columns=['gene',
                                  'function',
                                  'database',
                                  'db_ID',
                                  'gain',
                                  'loss',
                                  'NF1',
                                  'NF2',
                                  'NF3',
                                  'NF4',
                                  'NF5',
                                  'NF6',
                                  'NF7',
                                  'NF8',
                                  'NF9',
                                  'NF10',
                                  'NF11'])

for g_name, names in tqdm(group_idx.items()):
    sub_df = ori_df.loc[names, :]
    sub_g = sub_df.groupby('genename')
    result_df = result_df.append(pd.DataFrame.from_dict(
        {0: {'gene': g_name}}, orient='index'), sort=False)

    tmp_result_df = pd.DataFrame(columns=result_df.columns)
    for genename, ssub_g_idxs in sub_g.groups.items():
        ssub_df = sub_df.loc[ssub_g_idxs, :]
        row = ssub_df.iloc[0, :]
        new_dict = {'gene': genename.strip(),
                    'function': str(row['gene function description']).strip(),
                    'database': str(row['database']).strip() if isinstance(row['database'], str) else '',
                    'db_ID': str(row['KO number']).strip(),
                    'gain': ssub_df.loc[ssub_df.loc[:, 'gain or loss'] > 0, :].shape[0],
                    'loss': ssub_df.loc[ssub_df.loc[:, 'gain or loss'] < 0, :].shape[0],
                    'appear times': ssub_df.loc[ssub_df.loc[:, 'gain or loss'] != 0, :].shape[0]}
        for _, _row in ssub_df.iterrows():
            # _row = _row.astype(int)
            _cache = ','.join(
                list(map(lambda x: str(int(x)),
                         _row[['extant', 'ancestral', 'gain or loss']])))
            new_dict[_row['nod to nif']] = '' if _cache == '0,0,0' else _cache
        new_series = pd.DataFrame.from_dict({0: new_dict}, orient='index')

        tmp_result_df = tmp_result_df.append(new_series, sort=False)
    tmp_result_df = tmp_result_df.sort_values(['appear times', 'gain', 'loss'], ascending=False)
    result_df = result_df.append(tmp_result_df, sort=False)
    result_df = result_df.append(pd.DataFrame.from_dict(
        {0: {'gene': ''}}, orient='index'), sort=False)
    result_df = result_df.append(pd.DataFrame.from_dict(
        {0: {'gene': ''}}, orient='index'), sort=False)

result_df = result_df.drop('appear times', axis=1)
result_df.to_excel(output_xlsx, index=False)
