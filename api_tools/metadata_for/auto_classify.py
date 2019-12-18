import pandas as pd
from tqdm import tqdm
from os.path import join,expanduser
from api_tools.itol_func import *

fpath = expanduser("~/script/evolution_relative/api_tools/metadata_for/keyword.csv")


kw_df =  pd.read_csv(fpath,sep='\t',index_col=0)
kw_df.to_dict(orient='index')

def _classificated(ori_df):
    kw1='classification(auto)'
    kw2='habitat(auto)'
    kw3='matched keyword(auto)'

    for _,row in tqdm(ori_df.iterrows(),total=ori_df.shape[0]):
        row_text = ' ; '.join(map(str,row.values)).lower()
        _cache3 = set()
        _cache2 = set()
        for kw in keyword_mapping:
            l_kw = kw.lower()
            if l_kw in row_text:
                _cache3.add(kw)
                _cache2.add(keyword_mapping[kw])
            ori_df.loc[_,kw2] = ';'.join(sorted(_cache2))
            ori_df.loc[_,kw3] = ';'.join(sorted(_cache3))
        if 'metageno' in row_text:
            ori_df.loc[_,kw1] = 'MAGs'
        if 'single cell' in row_text:
            ori_df.loc[_,kw1] = 'SAGs'
        if 'Whole genome' in row_text or 'type strain' in row_text:
            ori_df.loc[_,kw1] = 'isolate'
    return ori_df

