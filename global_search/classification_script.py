import pandas as pd
from tqdm import tqdm

# classifiy_word_table = {'MAGs':['metageno'],
                        
#                         }

# infile = './nr_retrieve_hao/manual_annotated_biosample.csv'
def _classificated(ori_df):
#ori_df = pd.read_csv(infile,index_col=0)
    kw1='classification(auto)'
    kw2='habitat(auto)'
    ori_df.loc[:,'classification'] = ''
    ori_df.loc[:,'habitat'] = ''

    for _,row in tqdm(ori_df.iterrows(),total=ori_df.shape[0]):
        row_text = ';'.join(map(str,row.values)).lower()
        if 'metageno' in row_text:
            ori_df.loc[_,kw1] = 'MAGs'
        if 'single cell' in row_text:
            ori_df.loc[_,kw1] = 'SAGs'
        if 'soil' in row_text or 'terrestrial' in row_text or 'wetland' in row_text:
            ori_df.loc[_,kw2] = 'terrestrial'
        if 'sea' in row_text or 'marine' in row_text or 'ocean' in row_text:
            ori_df.loc[_,kw2] = 'marine'
        if 'wastewater' in row_text or 'activated sludge' in row_text:
            ori_df.loc[_,kw2] = 'waste water'
        if 'groundwater' in row_text:
            ori_df.loc[_,kw2] = 'ground water'
        if 'bioreactor' in row_text:
            ori_df.loc[_,kw2] = 'bioreactor'
        if 'bioreactor' in row_text and 'sludge' in row_text:
            ori_df.loc[_,kw2] = 'waste water'
        if 'freshwater' in row_text:
            ori_df.loc[_,kw2] = 'freshwater'
        if 'soda lake' in row_text:
            ori_df.loc[_,kw2] = 'alkaline aquatic'
        if 'subsurface' in row_text:
            ori_df.loc[_,kw2] = 'subsurface'
        if 'biofilter' in row_text or 'sand filter' in row_text:
            ori_df.loc[_,kw2] = 'artificial system'
        if 'Whole genome' in row_text or 'type strain' in row_text:
            ori_df.loc[_,kw1] = 'isolate'

        # if not pd.isna(row['attribute:host']) and str(row['attribute:host']) != 'not applicable':
        #     ori_df.loc[_,'habitat'] = 'host associated'
        if 'source' in row.index:
            if 'uncultured' in row['source'] or 'unidentified' in row['source']:
                ori_df.loc[_,kw2] = 'amplicons'
        list_c = list(ori_df.columns)
        list_c = list_c[:-2]
        list_c.insert(12,kw2)
        list_c.insert(12,kw1)
        ori_df = ori_df.reindex(columns=list_c)
    return ori_df
    # ori_df.to_csv('./nr_retrieve_hao/test.csv',index=1)