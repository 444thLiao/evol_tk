import pandas as pd
import sys
from os.path import exists,join
from tqdm import tqdm
from api_tools.itol_func import *  
from glob import glob 
from ete3 import Tree
import plotly.express as px


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
        

def modify_ID(now_dict,treeIDs):
    new_dict = {}
    failed_id = []
    for id in treeIDs:
        if id not in now_dict:
            if 'KU509367.1' in id:
                # manual mistake, it should be a nuc id
                _id = 'ANC58166.1'
            else:
                _id = [k for k in now_dict if k.split('.')[0] in id]
                if not _id:
                    failed_id.append(id)
                    continue
                else:
                    _id = _id[0]
            new_dict[id] = now_dict[_id]
        else:
            new_dict[id]= now_dict[id]
    return new_dict


def get_colors_general(ID2infos, now_info2style={}):
    colors = px.colors.qualitative.Dark24 + px.colors.qualitative.Light24
    remained_colors = [c for c in colors if c not in now_info2style.values()]
    info2style = {}
    for v in set(ID2infos.values()):
        one_color = remained_colors.pop(0)
        info2style.update({v: one_color})
    return ID2infos, info2style

outgroup_gene_names = {'K00370': ['dms', 'tor'],
                       'K00371': ['dms', ],
                       'K10535': ['nrfA', '_ONR'],
                       'K10944': ['bmo'],
                       'K10945': ['bmo'],
                       'K10946': ['bmo']}

ref_file = '/home-user/thliao/project/nitrogen_cycle/nitrification/reference_genomes/outgroup and reference.xlsx'
ref_df = pd.read_excel(ref_file, index_col=None)
ref_df = ref_df.loc[ref_df.loc[:, 'note'] != 'removed', :]
ref_id2habitat = {str(row['AA accession'])+'_'+str(row['gene name']): str(row['habitat (manual)'])
                      for _,row in ref_df.iterrows()}


gene_info = {'nxrA': 'K00370',
                      'nxrB': 'K00371',
                      'hao': 'K10535',
                      'amoA': 'K10944',
                      'amoB': 'K10945',
                      'amoC': 'K10946'}

if len(sys.argv) >= 2:
    file_list = sys.argv[1:]
    failed_f = []
    for fdir in tqdm(file_list):
        #try:
        gene = fdir.split('/')[0].split('_')[-1]
        ko = gene_info[gene]
        sub_ref_df = ref_df.loc[ref_df.loc[:, 'outgroup/ref for which KO'] == ko, :]
        sub_ref_df = sub_ref_df.loc[sub_ref_df.loc[:,
                                           'phylum/class'] != 'Thaumarchaeota', :]
        
        f = join(fdir,'full_info_new.xlsx')
        tfile = glob(join(fdir,'*.sorted.newick'))[0]
        tree = Tree(tfile,format=3)
        all_ids = list(tree.get_leaf_names())
        full_df = pd.read_excel(f,index_col=0)
        full_df = full_df.fillna('unknown')
        id2habitat = dict(zip(full_df.index,
                              full_df.loc[:,'habitat']))
        new_id2habitat = modify_ID(id2habitat,all_ids)
        new_id2habitat.update({k:v for k,v in ref_id2habitat.items() if k in all_ids})
        to_habitat_matrix(new_id2habitat,fdir)
        ####
        
        id2seq_type = {}
        for _,row in  full_df.iterrows():
            if row['BioSample'] == 'unknown':
                id2seq_type[_] = 'amplicons'
            else:
                id2seq_type[_] = 'with Genomes'
        seq_type_style = {'amplicons':'#0000ff',
                          'with Genomes':'#ff0000',}
        id2seq_type = modify_ID(id2seq_type,all_ids)
        matrix_text = to_matrix_shape(id2seq_type,'seq type',seq_type_style)
        with open(join(fdir,'seqtype_matrix.txt'),'w') as f1:
            f1.write(matrix_text)
        #### 
        id2info = {str(row['AA accession'])+'_'+str(row['gene name']): str(row['gene name'])
                      for _,row in sub_ref_df.iterrows()}
        
        id2info,info2col = get_colors_general(id2info,)
        
        to_color_labels_bg()
        write2colorbranch_clade(id2info,
                            dirname(tfile),
                            info2col,
                            treefile=tfile,
                            unique_id=ko,
                            info_name='branch_color',
                            no_legend=False)
        # if failed_id:
        #     print(fdir,failed_id)
    if failed_f:
        print(len(failed_f),' failed')
    else:
        print('all successfull')