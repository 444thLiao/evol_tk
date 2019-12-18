import pandas as pd
from tqdm import tqdm
from os.path import join

# classifiy_word_table = {'MAGs':['metageno'],
                        
#                         }

# infile = './nr_retrieve_hao/manual_annotated_biosample.csv'
def _classificated(ori_df):
#ori_df = pd.read_csv(infile,index_col=0)
    kw1='classification(auto)'
    kw2='habitat(auto)'
    ori_df.loc[:,kw1] = ''
    ori_df.loc[:,kw2] = ''

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
        if 'groundwater' in row_text:
            ori_df.loc[_,kw2] = 'ground water'
        if 'bioreactor' in row_text:
            ori_df.loc[_,kw2] = 'bioreactor'
        if ('bioreactor' in row_text and 'sludge' in row_text) or 'wastewater' in row_text or 'activated sludge' in row_text:
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
    
def diff_marine_non_marine(ori_df):
    ori_df = ori_df.copy()
    kw2='habitat(auto,diff marine/non-marine)'
    ori_df.loc[:,kw2] = ''

    for _,row in tqdm(ori_df.iterrows(),total=ori_df.shape[0]):
        row_text = ';'.join(map(str,row.values)).lower()
        if 'sea' in row_text or ('marine' in row_text and 'non-marine' not in row_text) or 'ocean' in row_text :
            ori_df.loc[_,kw2] = 'marine'
        elif 'soda lake' in row_text:
            ori_df.loc[_,kw2] = 'marine'
        elif 'lagoon' in row_text or 'brackish' in row_text:
            ori_df.loc[_,kw2] = 'marine'
        elif ('bioreactor' in row_text and 'sludge' in row_text) or 'wastewater' in row_text or 'activated sludge' in row_text:
            ori_df.loc[_,kw2] = 'marine'
        elif pd.isna(ori_df.loc[_,['habitat(auto)','habitat']]).all():
            ori_df.loc[_,kw2] = 'unknown'
        elif str(ori_df.loc[_,'habitat'])=='marine':
            ori_df.loc[_,kw2] = 'marine'
        elif str(ori_df.loc[_,'habitat'])=='non-marine':
            ori_df.loc[_,kw2] = 'non-marine'
        else:
            ori_df.loc[_,kw2] = 'non-marine'
        list_c = list(ori_df.columns)
        list_c = list_c[:-2]
        list_c.insert(12,kw2)
        #ori_df = ori_df.reindex(columns=list_c)
    return ori_df

if __name__ == "__main__":
    from api_tools.itol_func import *
    from ete3 import Tree
    metadata = './rawdata/plancto_final_habitat.xlsx'
    new_df = pd.read_excel(metadata)
    id2habitat = dict(zip(new_df.iloc[:,1] ,
                        new_df.loc[:,'habitat (marine-non-marine)']
                        )
                    )

    id2habitats = {k:[v] if not pd.isna(v) else ['unknown'] for k,v in id2habitat.items()}

    all_habitats = list(set([v[0] for k,v in id2habitats.items() ]))
    all_text = to_binary_shape(id2habitats,
                    {_:{'color':'#D68529'} if _=='non-marine' else {'color':'#0011FF'} for _ in all_habitats},
                    info_name='habitat',
                    omitted_other=True)
    with open(join('./itol_txt','general_habitat.txt'),'w') as f1:
        f1.write(all_text)
        
    # outgroup
    id2habitat.update({"GCA_000019665.1":'non-marine',
                       "GCA_000020225.1":'non-marine',
                       "GCA_000172155.1":'non-marine',
                       "GCA_001318295.1":'non-marine',
                       "GCA_001613545.1":'non-marine',
                       "GCA_900097105.1":'non-marine',
                       "GCA_001746835.1":'non-marine'})

    rows = []
    tree = Tree('../trees/iqtree/over20p_bac120.formatted.newick',format=3)
    for l in tree.get_leaves():
        name = l.name
        habitat = id2habitat.get(name,'-')
        if habitat == 'marine':
            rows.append('\t'.join([name,'1','0']))
            
        else:
            rows.append('\t'.join([name,'0','1']))
    with open('./m2nm.txt','w') as f1:
        f1.write('\n'.join(rows))

        
        