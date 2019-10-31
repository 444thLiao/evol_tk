import pandas as pd
import sys
from os.path import exists, join
from tqdm import tqdm
from api_tools.itol_func import *
from glob import glob
from ete3 import Tree
import plotly.express as px
from ete3 import NCBITaxa
ncbi = NCBITaxa()


def write2colorstrip(id2info, odir, info2color, unique_id, info_name='type',):
    content = to_color_strip(id2info, info2color, info_name=info_name)
    info_name = info_name.replace('/', '_')
    with open(join(odir, f'{unique_id}_{info_name}_colorstrip.txt'), 'w') as f1:
        f1.write(content)


def write2color_label_bg(id2info, odir, info2color, unique_id, info_name='type',):
    content = to_color_labels_bg(id2info, info2color, info_name=info_name)
    info_name = info_name.replace('/', '_')
    with open(join(odir, f'{unique_id}_{info_name}_color_label_bg.txt'), 'w') as f1:
        f1.write(content)


def write2colorbranch(id2info, odir, info2color, unique_id, info_name='type', no_legend=False):
    content = to_color_branch(
        id2info, info2color, dataset_name=info_name, no_legend=no_legend)
    info_name = info_name.replace('/', '_')
    with open(join(odir, f'{unique_id}_{info_name}_colorbranch.txt'), 'w') as f1:
        f1.write(content)


def write2colorbranch_clade(id2info, odir, info2color, treefile, unique_id, info_name='type', **kwargs):
    content = to_color_Clade(
        id2info, info2color, treefile, info_name, **kwargs)
    info_name = info_name.replace('/', '_')
    with open(join(odir, f'{unique_id}_{info_name}_colorbranch_clade.txt'), 'w') as f1:
        f1.write(content)


def write2binary_dataset(ID2infos, odir, info2style, unique_id):
    annotate_text = annotate_outgroup(ID2infos, info2style)
    with open(join(odir, f'{unique_id}_marker_outgroup_ref.txt'), 'w') as f1:
        f1.write(annotate_text)


habitat_colros = {'marine': '#0011FF',
                  'marine(symbiont)': '#0011FF',
                  'terrestrial': '#f99806',
                  'terrestrial(symbiont)': '#f99806',
                  'freshwater': '#88bbFF',
                  'waste water': '#50AF6E',
                  'artifical system':'#FA3705'
                  }


def to_habitat_matrix(id2habitat, fdir):
    matrix_text = to_matrix_shape(id2habitat, 'habitat', habitat_colros)
    with open(join(fdir, 'habitat_matrix.txt'), 'w') as f1:
        f1.write(matrix_text)


def modify_ID(now_dict, treeIDs):
    new_dict = {}
    failed_id = []
    for id in treeIDs:
        if id not in now_dict:
            if 'KU509367.1' in id:
                # manual mistake, it should be a nuc id
                _id = 'ANC58166.1'
                if _id not in now_dict:
                    continue
            else:
                _id = [k for k in now_dict if k.split('.')[0] in id]
                if not _id:
                    failed_id.append(id)
                    continue
                else:
                    _id = _id[0]
            new_dict[id] = now_dict[_id]
        else:
            new_dict[id] = now_dict[id]
    return new_dict


def get_colors_general(ID2infos, now_info2style={}):
    colors = px.colors.qualitative.Dark24 + px.colors.qualitative.Light24
    remained_colors = [c for c in colors if c not in now_info2style.values()]
    info2style = {}
    for v in set(ID2infos.values()):
        if v in now_info2style:
            info2style.update({v: now_info2style[v]})

        else:
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
                  for _, row in ref_df.iterrows()}


gene_info = {'nxrA': 'K00370',
             'nxrB': 'K00371',
             'hao': 'K10535',
             'amoA': 'K10944',
             'amoB': 'K10945',
             'amoC': 'K10946'}

file_list = ['nr_retrieve_nxrA/cluster_95_filtered_lengths.fa_aln.dir/iqtree.treefile',
             'nr_retrieve_hao/filtered_by_kegg.faa_aln.dir/iqtree.no_trim.treefile/',
             'nr_retrieve_amoB/filtered_by_kegg.faa_aln.dir/iqtree.treefile',
             'nr_retrieve_amoC/filtered_by_kegg.faa_aln.dir/iqtree.treefile',
             'nr_retrieve_removeENV_amoA/cluster_98_aln.dir/iqtree.treefile']
# './nr_retrieve_removeENV_true_amoA/removed_ENV.faa_aln.dir/iqtree.treefile' 
if len(sys.argv) >= 2:
    file_list = sys.argv[1:]
    failed_f = []
    for fdir in tqdm(file_list):
        # try:
        gene = fdir.split('/')[0].split('_')[-1]
        ko = gene_info[gene]
        sub_ref_df = ref_df.loc[ref_df.loc[:,
                                           'outgroup/ref for which KO'] == ko, :]
        sub_ref_df = sub_ref_df.loc[sub_ref_df.loc[:,
                                                   'phylum/class'] != 'Thaumarchaeota', :]

        f = join(fdir, 'full_info_new.xlsx')
        if not exists(f):
            f = join(fdir, 'full_info.xlsx')
        tfile = glob(join(fdir, '*.sorted.newick'))[0]
        tree = Tree(tfile, format=3)
        all_ids = list(tree.get_leaf_names())
        full_df = pd.read_excel(f, index_col=0)
        full_df = full_df.fillna('unknown')
        id2habitat = dict(zip(full_df.index,
                              full_df.loc[:, 'habitat']))
        new_id2habitat = modify_ID(id2habitat, all_ids)
        new_id2habitat.update(
            {k: v for k, v in ref_id2habitat.items() if k in all_ids})
        to_habitat_matrix(new_id2habitat, fdir)
        #### seq source

        id2seq_type = {}
        for _, row in full_df.iterrows():
            if row['BioSample'] == 'unknown':
                id2seq_type[_] = 'amplicons'
            else:
                id2seq_type[_] = 'with Genomes'
        seq_type_style = {'amplicons': '#0000ff',
                          'with Genomes': '#ff0000', }
        id2seq_type = modify_ID(id2seq_type, all_ids)
        content = to_color_strip(id2seq_type, seq_type_style,'seq type')
        with open(join(fdir, 'seqtype_colorstrip.txt'), 'w') as f1:
            f1.write(content)
        ######## gene name 
        
        gname = ['bmo','pmo','pxm','amo']
        gcolors = '#2E91E5,#E15F99,#1CA71C,#FB0D0D'.split(',')
        
        id2info = {str(row['AA accession'])+'_'+str(row['gene name']): str(row['gene name'])
                   for _, row in sub_ref_df.iterrows()}
        infos = set(id2info.values())
        now_colors_dict = {}
        for name in infos:
            if name[:-1] in gname:
                now_colors_dict[name] = dict(zip(gname,gcolors))[name[:-1]]
        
        id2info, info2col = get_colors_general(id2info,now_info2style=now_colors_dict)
        id2info = modify_ID(id2info, all_ids)
        template_text = to_color_strip(
            id2info, info2col, info_name='gene name')
        with open(join(fdir, 'gene_name_strip.txt'), 'w') as f1:
            f1.write(template_text)
        ############## taxonomy
        id2tax = {}
        for aid, tid in zip(full_df.index, full_df.loc[:, 'taxid']):
            if tid == 'unknown':
                id2tax[aid] = 'unknown'
                continue
            lineage = ncbi.get_lineage(tid)
            rank = ncbi.get_rank(lineage)
            rank = {v: k for k, v in rank.items()}
            names = ncbi.get_taxid_translator(lineage)
            if names.get(rank.get('phylum', ''), 'ENV') == 'Proteobacteria':
                id2tax[aid] = names.get(rank.get('class', ''), 'ENV')
            else:
                id2tax[aid] = names.get(rank.get('phylum', ''), 'ENV')

        id2tax = {k: 'CPR' if 'candidat' in v.lower() else v
                  for k, v in id2tax.items()}
        Interested_tax = {'Thaumarchaeota': '#358f0f',
                          'Nitrospirae': '#edc31d',
                          'Chloroflexi': '#e41a1c',
                          'Gammaproteobacteria': '#7ffcfa',
                          'Betaproteobacteria': '#956cb4',
                          'Alphaproteobacteria': '#8c613c',
                          'Actinobacteria': '#11FF11',
                          'Planctomycetes': '#ff44bb',
                          'ENV': '#B54B4A',
                          'CPR': '#74A45B'
                          }
        id2tax = {k: v for k, v in id2tax.items() if v in Interested_tax}
        id2tax = modify_ID(id2tax, all_ids)
        id2info, info2col = get_colors_general(
            id2tax, now_info2style=Interested_tax)
        # write2colorbranch_clade(id2info,
        #                         dirname(tfile),
        #                         info2col,
        #                         treefile=tfile,
        #                         unique_id=ko,
        #                         info_name='branch_color',
        #                         no_legend=False)
        template_text = to_color_strip(
            id2info, info2col, info_name='gene name')
        with open(join(fdir, 'tax_strip.txt'), 'w') as f1:
            f1.write(template_text)

        # full_text = to_color_labels_bg(id2info,info2col,info_name='part tax')
        # with open(join(fdir, 'bg_phylum.txt'),'w') as f1:
        #     f1.write(full_text)
        # if failed_id:
        #     print(fdir,failed_id)
    if failed_f:
        print(len(failed_f), ' failed')
    else:
        print('all successfull')
