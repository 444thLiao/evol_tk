import pandas as pd
import sys
from os.path import exists,join
from tqdm import tqdm
from api_tools.itol_func import *  
from glob import glob 
from ete3 import Tree
def write2colorstrip(id2info,odir,info2color, unique_id,info_name='type',):
    content = to_color_strip(id2info,info2color,info_name=info_name)
    info_name = info_name.replace('/','_')
    with open(join(odir, f'{unique_id}_{info_name}_colorstrip.txt'), 'w') as f1:
        f1.write(content)

def write2color_label_bg(id2info,odir,info2color, unique_id,info_name='type',):
    content = to_color_labels_bg(id2info,info2color,info_name=info_name)
    info_name = info_name.replace('/','_')
    with open(join(odir, f'{unique_id}_{info_name}_color_label_bg.txt'), 'w') as f1:
        f1.write(content)
        
def write2colorbranch(id2info,odir,info2color, unique_id,info_name='type',no_legend=False):
    content = to_color_branch(id2info,info2color,dataset_name=info_name,no_legend=no_legend)
    info_name = info_name.replace('/','_')
    with open(join(odir, f'{unique_id}_{info_name}_colorbranch.txt'), 'w') as f1:
        f1.write(content)

def write2colorbranch_clade(id2info,odir,info2color,treefile, unique_id,info_name='type',**kwargs):
    content = to_color_Clade(id2info,info2color,treefile,info_name,**kwargs)
    info_name = info_name.replace('/','_')
    with open(join(odir, f'{unique_id}_{info_name}_colorbranch_clade.txt'), 'w') as f1:
        f1.write(content)


def write2binary_dataset(ID2infos, odir,info2style, unique_id):
    annotate_text = annotate_outgroup(ID2infos,info2style)
    with open(join(odir,f'{unique_id}_marker_outgroup_ref.txt'),'w') as f1:
        f1.write(annotate_text)


habitat_colros = {'marine':'#0099FF',
                  'marine(symbiont)':'#0099FF',
                  'terrestrial':'#f99806',
                  'terrestrial(symbiont)':'#f99806',
                  'waste water':'#aaffaa',
                  }   
def to_habitat_matrix(id2habitat,fdir):
    matrix_text = to_matrix_shape(id2habitat,'habitat',habitat_colros)
    with open(join(fdir,'habitat_matrix.txt'),'w') as f1:
        f1.write(matrix_text)
        
        
if len(sys.argv) >= 2:
    file_list = sys.argv[1:]
    failed_f = []
    for fdir in tqdm(file_list):
        #try:
        f = join(fdir,'full_info_new.xlsx')
        t = glob(join(fdir,'*.sorted.newick'))[0]
        tree = Tree(t,format=3)
        all_ids = list(tree.get_leaf_names())
        full_df = pd.read_excel(f,index_col=0)
        full_df = full_df.fillna('unknown')
        id2habitat = dict(zip(full_df.index,
                              full_df.loc[:,'habitat']))
        new_id2habitat = {}
        failed_id = []
        for id in all_ids:
            if id not in id2habitat:
                if 'KU509367.1' in id:
                    # manual mistake, it should be a nuc id
                    _id = 'ANC58166.1'
                else:
                    _id = [k for k in id2habitat if k.split('.')[0] in id]
                    if not _id:
                        failed_id.append(id)
                        continue
                    else:
                        _id = _id[0]
                new_id2habitat[id] = id2habitat[_id]
            else:
                new_id2habitat[id]= id2habitat[id]
        to_habitat_matrix(new_id2habitat,fdir)
        ####
        
        id2seq_type = {}
        for _,row in  full_df.iterrows():
            if row['BioSample'] == 'unknown':
                id2seq_type[_] = 'amplicons'
            
        
        if failed_id:
            print(fdir,failed_id)
    if failed_f:
        print(len(failed_f),' failed')
    else:
        print('all successfull')