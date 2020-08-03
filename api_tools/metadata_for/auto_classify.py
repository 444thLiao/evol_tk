import pandas as pd
from tqdm import tqdm
from os.path import join,expanduser
# from api_tools.itol_func import *

fpath = expanduser("~/script/evolution_relative/api_tools/metadata_for/keyword.csv")


kw_df = pd.read_csv(fpath,sep='\t',index_col=0)
kw_dict = kw_df.to_dict(orient='index')

def _classificated(ori_df):
    kw1='classification(auto)'
    kw3='matched keyword(auto)'
    kws1='habitat_set1(auto)'
    kws2='habitat_set2(auto)'
    kws3='habitat_set3(auto)'
    subset_columns = [_ for _ in ori_df.columns if _ not in {kw1,kw3,kws1,kws2,kws3}]
    ori_df = ori_df.loc[:,subset_columns]
    for _,row in tqdm(ori_df.iterrows(),total=ori_df.shape[0]):
        row_text = ' ; '.join(map(str,row.values)).lower()
        _cache3 = set()
        _caches1 = set()
        _caches2 = set()
        _caches3 = set()
        for kw,sub_d in kw_dict.items():
            l_kw = kw.lower()
            if l_kw in row_text:
                _cache3.add(kw)
                _caches1.add(sub_d['set1']) if not pd.isna(sub_d['set1']) else None
                _caches2.add(sub_d['set2']) if not pd.isna(sub_d['set2']) else None
                _caches3.add(sub_d['set3']) if not pd.isna(sub_d['set3']) else None
            ori_df.loc[_,kws1] = ';'.join(sorted(_caches1))
            ori_df.loc[_,kws2] = ';'.join(sorted(_caches2))
            ori_df.loc[_,kws3] = ';'.join(sorted(_caches3))
            ori_df.loc[_,kw3] = ';'.join(sorted(_cache3))
        if 'metageno' in row_text:
            ori_df.loc[_,kw1] = 'MAGs'
        if 'single cell' in row_text:
            ori_df.loc[_,kw1] = 'SAGs'
        if 'whole genome' in row_text or 'type strain' in row_text or 'whole-genome' in row_text or 'isolate' in row_text:
            ori_df.loc[_,kw1] = 'isolate'
        ori_df.loc[:,kw1].fillna('unidentified',inplace=True)
    return ori_df

