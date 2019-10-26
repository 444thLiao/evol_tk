import pandas as pd
source_file = '/home-user/thliao/resource/NCBI2habitat.csv'
import sys
from os.path import exists,join

if len(sys.argv) >= 2:
    file_list = sys.argv[1:]
    source_df = pd.read_csv(source_file, encoding='GBK')
    source_df.loc[:, 'tmp'] = source_df.iloc[:, 0] + \
        ';'+source_df.iloc[:, 1]
    source_df = source_df.set_index('tmp')
    source_df = source_df.loc[~source_df.index.duplicated(), :]
    
    for fdir in file_list:
        f = join(fdir,'full_info.xlsx')
        full_df = pd.read_excel(f,index_col=0)
        if exists(source_file):

            # bioproject
            tmp = [';'.join(list(map(str, row))[:2])
                    for row in full_df.loc[:, ['BioProject', 'BioSample']].values]

            _d1 = source_df.reindex(tmp)
            for idx, (_, v) in enumerate(full_df.iterrows()):
                if not pd.isna(_d1.iloc[idx, 2]) and pd.isna(v['habitat']):
                    # print(_,_d1.iloc[idx,2])
                    full_df.loc[_, 'habitat'] = _d1.iloc[idx, 2]
        full_df.to_excel(f.replace('.xlsx','_new.xlsx'),
                            index=1, index_label='protein accession')