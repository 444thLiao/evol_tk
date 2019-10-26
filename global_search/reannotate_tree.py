import pandas as pd
import sys
from os.path import exists,join
from tqdm import tqdm
from api_tools.itol_func import *   
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
                  'terrestrial':'#f99806',
                  'waste water':'#aaffaa',
                  }   
def to_habitat_matrix(id2habitat,fdir):
    all_habitat = id2habitat
    matrix_text = to_matrix_shape(id2habitat,'habitat',habitat_colros)
    with open(join(fdir,'habitat_matrix.txt'),'w') as f1:
        f1.write(matrix_text)
        
        
if len(sys.argv) >= 2:
    file_list = sys.argv[1:]
    failed_f = []
    for fdir in tqdm(file_list):
        try:
            f = join(fdir,'full_info_new.xlsx')
            full_df = pd.read_excel(f,index_col=0)
            full_df = full_df.fillna('Unknown')
            id2habitat = dict(zip(full_df.index,full_df.loc[:,'habitat']))
            to_habitat_matrix(id2habitat,fdir)
        except:
            failed_f.append(fdir)
    if failed_f:
        print(len(failed_f),' failed')
    else:
        print('all successfull')