import pandas as pd
from tqdm import tqdm
from os.path import join,expanduser
from api_tools.itol_func import *

fpath = expanduser("~/script/evolution_relative/api_tools/metadata_for/keyword.csv")
keyword_mapping =  dict([row.strip('\n').splti('\t') for row in open(fpath).readlines() if row])

def _classificated(ori_df):
    # kw1='classification(auto)'
    kw2='habitat(auto)'
    kw3='collected keyword(auto)'
    ori_df.loc[:,kw3] = ''
    ori_df.loc[:,kw2] = ''

    for _,row in tqdm(ori_df.iterrows(),total=ori_df.shape[0]):
        row_text = ' ; '.join(map(str,row.values)).lower()
        _cache3 = set()
        _cache2 = set()
        for kw in keyword_mapping:
            if kw in row_text:
                _cache3.add(kw)
                _cache2.add(keyword_mapping[kw])
            ori_df.loc[_,'kw2'] = ';'.join(sorted(_cache2))
            ori_df.loc[_,'kw3'] = ';'.join(sorted(_cache3))
        # if 'metageno' in row_text:
        #     ori_df.loc[_,kw1] = 'MAGs'
        # if 'single cell' in row_text:
        #     ori_df.loc[_,kw1] = 'SAGs'
        if 'Whole genome' in row_text or 'type strain' in row_text:
            ori_df.loc[_,kw1] = 'isolate'

        # if not pd.isna(row['attribute:host']) and str(row['attribute:host']) != 'not applicable':
        #     ori_df.loc[_,'habitat'] = 'host associated'
        
    return ori_df
    
def diff_marine_non_marine(ori_df):
    ori_df = ori_df.copy()
    kw2='habitat(auto,diff marine/non-marine)'
    ori_df.loc[:,kw2] = ''

    for _,row in tqdm(ori_df.iterrows(),
                      total=ori_df.shape[0]):
        row = row[[_ for _ in row.index if not '(auto)' in _]]
        row_text = ';'.join(map(str,row.values)).lower()
        if ('marine' in row_text and 'non-marine' not in row_text) or 'ocean' in row_text :
            ori_df.loc[_,kw2] = 'marine'
        elif 'soda lake' in row_text:
            ori_df.loc[_,kw2] = 'marine'
        elif 'lagoon' in row_text or 'brackish' in row_text:
            ori_df.loc[_,kw2] = 'marine'
        elif ('bioreactor' in row_text and 'sludge' in row_text) or 'wastewater' in row_text or 'activated sludge' in row_text:
            ori_df.loc[_,kw2] = 'marine'
        elif pd.isna(ori_df.loc[_,['habitat(auto)']]).all():
            ori_df.loc[_,kw2] = 'unknown'
        # elif str(ori_df.loc[_,'habitat'])=='marine':
        #     ori_df.loc[_,kw2] = 'marine'
        # elif str(ori_df.loc[_,'habitat'])=='non-marine':
        #     ori_df.loc[_,kw2] = 'non-marine'
        else:
            ori_df.loc[_,kw2] = 'non-marine'
        list_c = list(ori_df.columns)
        list_c = list_c[:-2]
        list_c.insert(12,kw2)
        #ori_df = ori_df.reindex(columns=list_c)
    return ori_df
