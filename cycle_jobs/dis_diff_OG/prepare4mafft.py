import json
import os
from os.path import join, exists,basename,dirname
from glob import glob
import pandas as pd
from subprocess import check_call
from Bio import SeqIO
import seaborn as sns
from ete3 import Tree
import plotly.express as px

# load all necessary data
odir = 'json_dump_v2'
with open(join(odir, 'ko2og.json'), 'r') as f1:
    ko2og = json.load(f1)
with open(join(odir, 'manually_blast_r.json'), 'r') as f1:
    manually_blast_r = json.load(f1)
with open(join(odir, 'ko2og2names.json'), 'r') as f1:
    ko2og2names = json.load(f1)

indir = 'genome_protein_files_more'
og_tsv = f'./{indir}/OrthoFinder/Results_Oct01/Orthogroups/Orthogroups.tsv'
og_df = pd.read_csv(og_tsv, sep='\t', low_memory=False, index_col=0)
def rename(x):
    if pd.isna(x):
        return x
    if ', ' not in x:
        return str(x).split(' ')[0]
    else:
        return ', '.join([str(_).split(' ')[0] for _ in x.split(', ')])
og_df = og_df.applymap(rename)

genome_info = './genome_info_full.xlsx'
g_df = pd.read_excel(genome_info, index_col=0)

# name2dirname
name2dirname = {}
for ofile in glob(join(indir,'*.faa')):
    real_ofile = os.path.realpath(ofile)
    name = basename(ofile).replace('.faa','')
    if '/prokka_o/' in real_ofile:
        dir_name = basename(dirname(real_ofile)).replace('.faa','')
        name2dirname[name] = dir_name

# special_annotate_file (out group )
# generate annoating file
def generate_id2org(og_names, ofile):
    all_seq_ids = [_.id for _ in SeqIO.parse(ofile, format='fasta')]
    if isinstance(og_names,str):
        og_names = [og_names]
    id2org = {}
    for og_name in og_names:
        row = og_df.loc[og_name, :]
        for seq_id in all_seq_ids:
            org_ids = row.index[row.fillna('').str.contains(seq_id)]
            if not list(org_ids):
                #print('Error', seq_id, og_name)
                continue
            org_id = org_ids[0]
            id2org[seq_id] = org_id
    return id2org

def get_color_info(each_og, ofile,info_col='type',extra={}):
    id2org = generate_id2org(each_og, ofile)
    colors = px.colors.qualitative.Dark24 + px.colors.qualitative.Light24
    id2info = {}
    for id, org in id2org.items():
        org = name2dirname.get(org,org)
        if org in g_df.index:
            #name = g_df.loc[org, 'genome name']
            id2info[id] = g_df.loc[org,info_col]
        elif info_col == 'phylum/class' and name2dirname.get(org,org) in sname2info:
            id2info[id] = sname2info[org]
    if extra:
        id2info.update(extra)   
    set_v = set(id2info.values())
    num_v = len(set_v)
    cols = colors[:num_v]
    total_info2col = dict(zip(set_v,cols))
    if info_col in color_scheme:
        _info2col = color_scheme[info_col]
        info2col = {k:_info2col[k] for k in set_v if k in _info2col}
        fix_missed = {k:total_info2col[k] for k in set_v if k not in info2col}
        info2col.update(fix_missed)
    else:
        info2col = total_info2col.copy()
    return id2info,info2col

# rename
def to_label(each_og,ofile,ko_name=None):
    template_text = open(
        '/home-user/thliao/template_txt/labels_template.txt').read()
    id2org = generate_id2org(each_og, ofile)
    full_text = template_text[::]
    for id, org in id2org.items():
        if org in g_df.index:
            name = g_df.loc[org, 'genome name']
        else:
            name = id
        full_text += '%s,%s\n' % (id, name)
    if ko_name is not None:
        each_og = ko_name
    with open(join(odir, f'label_{each_og}.txt'), 'w') as f1:
        f1.write(full_text)

# rename the internal node
def renamed_tree(in_tree_file,outfile):
    count = 0
    t = Tree(open(in_tree_file).read())
    for n in t.traverse():
        if not n.name:
            n.name = 'Internal_%s' % count
            count += 1
    t.write(outfile=outfile,format=3)
    return t

color_scheme = {'type':{'NOB': '#e41a1c', 'comammox': '#edc31d', 
                        'AOB': '#bad5b9', 'AOA': '#358f0f'},
                'phylum/class':{'Thaumarchaeota': '#358f0f',
                                'Nitrospirae': '#edc31d',
                                'Gammaproteobacteria': '#78fce0',
                                'Chloroflexi': '#e41a1c',
                                'Betaproteobacteria': '#956cb4',
                                'Alphaproteobacteria': '#8c613c'}

                }

# annotate MAGs lineage
remained_ID = og_df.columns.difference(g_df.index)
MAG_annotate_file = '/home-user/thliao/data/metagenomes/update_0928_nitrification/confirmed_locus2info.tsv'
MAG_annotate_df = pd.read_csv(MAG_annotate_file,sep='\t',index_col=0)
derep_df = MAG_annotate_df.drop_duplicates('sample name')
filtered_df = derep_df.loc[derep_df.loc[:,'sample name'].isin(remained_ID),:]
phylum_from_metadata_count = filtered_df.groupby('phylum(from metadata)').count().iloc[:,0]
sname2phylum_metadata = dict(zip(MAG_annotate_df.loc[:,'sample name'],MAG_annotate_df.loc[:,'phylum(from metadata)']))
sname2class_metadata = dict(zip(MAG_annotate_df.loc[:,'sample name'],MAG_annotate_df.loc[:,'class(from metadata)']))

sname2phylum_metadata = {k:v for k,v in sname2phylum_metadata.items() if not pd.isna(v)}
sname2info = {k:v if v !='Proteobacteria' else sname2class_metadata[k] for k,v in sname2phylum_metadata.items()}
sname2info = {k:v for k,v in sname2info.items() if not pd.isna(v)}


# data dependent transform

def deduced_legend(info2color,info_name='dataset',sep=','):
    
    legend_title = info_name
    legend_shape = sep.join(['1'] * len(info2color))
    legend_colors = sep.join([_
                              for _ in info2color.values()])
    legend_labels = sep.join(list(info2color.keys()))
    legend_text = f"""
LEGEND_TITLE{sep}{legend_title}
LEGEND_SHAPES{sep}{legend_shape}
LEGEND_COLORS{sep}{legend_colors}
LEGEND_LABELS{sep}{legend_labels}"""
    return legend_text

def to_color_strip(id2info,info2color,info_name='dataset'):
    template_text = open(
        '/home-user/thliao/template_txt/dataset_color_strip_template.txt').read()
    id2col = {id:info2color[info] for id,info in id2info.items()}
    annotate_text = '\n'.join(['%s,%s\n' % (id,col) for id,col in id2col.items()])
    legend_text = deduced_legend(info2color,info_name)
    
    template_text = template_text.format(legend_text=legend_text,
                     dataset_label=info_name)
    info_name = info_name.replace('/','_')
    return template_text+'\n'+annotate_text
    
def to_color_branch(ID2info,info2color,dataset_name='color branch'):
    # clade for 
    template_text = open(
        '/home-user/thliao/template_txt/dataset_styles_template.txt').read()
    id2col = {ID:info2color[info] for ID,info in ID2info.items()}
    each_template = '{ID}\t{TYPE}\t{WHAT}\t{COLOR}\t{WIDTH_OR_SIZE_FACTOR}\t{STYLE}\t{BACKGROUND_COLOR}\n'
    legend_text = deduced_legend(info2color,dataset_name,sep='\t')
    
    template_text = template_text.format(dataset_label=dataset_name,
                                         legend_text=legend_text)
    rows = [each_template.format(ID=ID,
                                 TYPE='branch',
                                 WHAT='node',
                                 COLOR=color,
                                 WIDTH_OR_SIZE_FACTOR=7,
                                 STYLE='normal',
                                 BACKGROUND_COLOR='')
        for ID,color in id2col.items()]
    rows += [each_template.format(ID=ID,
                                 TYPE='label',
                                 WHAT='node',
                                 COLOR=color,
                                 WIDTH_OR_SIZE_FACTOR=1,
                                 STYLE='bold',
                                 BACKGROUND_COLOR='')
        for ID,color in id2col.items()]
    return template_text + '\n'.join(rows)

def to_color_Clade(ID2info,info2color,tree,dataset_name='color branch clade'):
    # clade for 
    template_text = open(
        '/home-user/thliao/template_txt/dataset_styles_template.txt').read()
    id2col = {ID:info2color[info] for ID,info in ID2info.items()}
    each_template = '{ID}\t{TYPE}\t{WHAT}\t{COLOR}\t{WIDTH_OR_SIZE_FACTOR}\t{STYLE}\t{BACKGROUND_COLOR}\n'
    legend_text = deduced_legend(info2color,dataset_name,sep='\t')
    
    tree_obj = Tree(tree,format=3)
    def collapsed_leaf(node):
        if node.is_leaf():
            return True
        else:
            leafs_all = node.get_leaves()
            children_names = [ID2info.get(_.name,'') for _ in leafs_all]
            if '' in children_names:
                return False
            if len(set(children_names)) == 1:
                return True
            else:
                return False
        return False
    new_tree_obj = Tree(tree_obj.write(is_leaf_fn=collapsed_leaf))
    new_leaves_names = [_.name for _ in new_tree_obj.get_leaves()]
    internal_nodes = [tree_obj.search_nodes(name =_)[0] for _ in new_leaves_names]
    internal_nodes = [_ for _ in internal_nodes if not _.is_leaf()]
    internal_node2info = {n.name:info2color.get(ID2info.get(n.get_leaves()[0].name,''),'') 
                          for n in internal_nodes}
    internal_node2info = {n.name:info2color[ID2info[n.get_leaves()[0].name]] 
                          for n in internal_nodes}

    template_text = template_text.format(dataset_label=dataset_name,
                                         legend_text=legend_text)
    rows = [each_template.format(ID=ID,
                                 TYPE='branch',
                                 WHAT='node',
                                 COLOR=color,
                                 WIDTH_OR_SIZE_FACTOR=7,
                                 STYLE='normal',
                                 BACKGROUND_COLOR='')
        for ID,color in id2col.items()]
    rows += [each_template.format(ID=ID,
                                 TYPE='label',
                                 WHAT='node',
                                 COLOR=color,
                                 WIDTH_OR_SIZE_FACTOR=1,
                                 STYLE='bold',
                                 BACKGROUND_COLOR='')
        for ID,color in id2col.items()]
    rows += [each_template.format(ID=ID,
                                 TYPE='branch',
                                 WHAT='clade',
                                 COLOR=color,
                                 WIDTH_OR_SIZE_FACTOR=7,
                                 STYLE='normal',
                                 BACKGROUND_COLOR='')
        for ID,color in internal_node2info.items()]
    return template_text + '\n'.join(rows)        
    
def write2colorstrip(id2info,info2color, unique_id,info_name='type',):
    content = to_color_strip(id2info,info2color,info_name=info_name)
    info_name = info_name.replace('/','_')
    with open(join(odir, f'{unique_id}_{info_name}_colorstrip.txt'), 'w') as f1:
        f1.write(content)

def write2colorbranch(id2info,info2color, unique_id,info_name='type',):
    content = to_color_branch(id2info,info2color,dataset_name=info_name)
    info_name = info_name.replace('/','_')
    with open(join(odir, f'{unique_id}_{info_name}_colorbranch.txt'), 'w') as f1:
        f1.write(content)

def write2colorbranch_clade(id2info,info2color,treefile, unique_id,info_name='type',):
    content = to_color_Clade(id2info,info2color,treefile,info_name)
    info_name = info_name.replace('/','_')
    with open(join(odir, f'{unique_id}_{info_name}_colorbranch_clade.txt'), 'w') as f1:
        f1.write(content)
       
# ref or outgroup seq, additionally add to
ref_file = '/home-user/thliao/project/nitrogen_cycle/nitrification/reference_genomes/outgroup and reference.xlsx'
ref_df = pd.read_excel(ref_file,index_col=None)
ref_df = ref_df.loc[ref_df.loc[:,'note']!='removed',:] 
def get_add_text(sub_df,used_ids):
    new_ref = []
    t_text = '' 
    id2info = {}
    for _,row in sub_df.iterrows():
        aa_id = row['AA accession']
        gene_name = row['gene name']
        seq = row['seq']
        info = row['phylum/class']
        if aa_id not in used_ids:
            t_text+= f'>{aa_id}_{gene_name}\n{seq}\n'
            
            id2info[f'{aa_id}_{gene_name}'] = info
        else:
            new_ref.append(aa_id)
    return t_text,new_ref,id2info

def annotate_outgroup(sub_df,ref_others=[]):
    template_text = open(
        '/home-user/thliao/template_txt/dataset_binary_template.txt').read()
    template_text = template_text.format(filed_shape='3',
                                         filed_label='outgroup/reference')
    annotate_text = ''
    for _,row in sub_df.iterrows():
        aa_id = row['AA accession']
        gene_name = row['gene name']
        name = f'{aa_id}_{gene_name}'
        if row['type'] == 'outgroup':
            annotate_text += '\t'.join([name,'0']) +'\n'
        else:
            annotate_text += '\t'.join([name,'1']) +'\n'
    for rid in ref_others:
        annotate_text += '\t'.join([rid,'1']) +'\n'
    return template_text + annotate_text


odir = join('./align_v2', 'complete_ko')
os.makedirs(odir, exist_ok=1)
for ko,og_list in ko2og.items():
    sub_ref_df = ref_df.loc[ref_df.loc[:,'outgroup/ref for which KO']==ko,:]
    predir = dirname(dirname(og_tsv))
    fa_files = [f'{predir}/Orthogroup_Sequences/{each_og}.fa' for each_og in og_list]
    used_ids = [record.id for fa_file in fa_files for record in SeqIO.parse(fa_file,format='fasta')]
    
    add_text,used_ref_ids,ref_id2info = get_add_text(sub_ref_df,used_ids)
    annotate_text = annotate_outgroup(sub_ref_df,ref_others=used_ref_ids)
    new_file = join(odir, ko+'.fa')
    if not exists(new_file):
        fa_file = ' '.join(fa_files)
        check_call(
            f'cat {fa_file} > {new_file}', shell=1)
    ori_text = open(new_file,'r').read()
    with open(new_file,'w') as f1:
        f1.write(add_text+ori_text)
    
    ofile = join(odir, ko+'.aln')
    if not exists(ofile):
        check_call(
            f'mafft --anysymbol --thread -1 {new_file} > {ofile}', shell=1)
    if not exists(ofile.replace('.aln','.treefile')):
       check_call(f'iqtree -nt 64 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -bb 1000 -pre {odir}/{ko} -s {ofile}',shell=1)
    #renamed_tree(ofile.replace('.aln','.treefile'),
    #             ofile.replace('.aln','.newick'))
    
    to_label(og_list,ofile,ko_name=ko)
    # annotate with type
    id2info,info2col = get_color_info(og_list,ofile,info_col='type')
    write2colorstrip(id2info,info2col,unique_id=ko,info_name='type')
    # annotate with phylum/class as a color strip
    id2info,info2col = get_color_info(og_list,ofile,info_col='phylum/class',extra=ref_id2info)
    write2colorstrip(id2info,info2col,unique_id=ko,info_name='phylum/class')
    # annotate with tree
    #write2colorbranch(id2info,info2col,unique_id=ko,info_name='branch_color')
    write2colorbranch_clade(id2info,
                            info2col,
                            treefile=ofile.replace('.aln','.newick'),
                            unique_id=ko,info_name='branch_color')
    # annotate the marker and outgroup
    with open(join(odir,f'{ko}_marker_outgroup_ref.txt'),'w') as f1:
        f1.write(annotate_text)
        

# # separated OG align
# odir = join('./align', 'single_OG')
# os.makedirs(odir, exist_ok=1)
# for ko, og_list in ko2og.items():
#     sub_ref_df = ref_df.loc[ref_df.loc[:,'outgroup/ref for which KO']==ko,:]
#     for each_og in og_list:
#         fa_file = f'./genome_protein_files/OrthoFinder/Results_Sep27/Orthogroup_Sequences/{each_og}.fa'
#         used_ids = [record.id for record in SeqIO.parse(fa_file,format='fasta')]
#         add_text,used_ref_ids = get_add_text(sub_ref_df,used_ids)
#         annotate_text = annotate_outgroup(sub_ref_df,ref_others=used_ref_ids)
#         new_fa_file = join(odir, each_og+'.fa')
#         with open(new_fa_file,'w') as f1:
#             f1.write(add_text+open(fa_file,'r').read())
#         ofile = join(odir, each_og+'.aln')
#         if not exists(ofile):
#             check_call(
#                 f'mafft --anysymbol --thread -1 {new_fa_file} > {ofile}', shell=1)
#         if not exists(ofile.replace('.aln','.treefile')):
#             check_call(f'iqtree -nt 32 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre {odir}/{each_og} -s {ofile}',shell=1)
#         to_label(each_og,ofile)
#         to_color_strip(each_og,ofile,info_col='type')
#         to_color_strip(each_og,ofile,info_col='phylum/class')
#         with open(join(odir,f'marker_{each_og}_outgroup_ref.txt'),'w') as f1:
#             f1.write(annotate_text)


# special process for differential AOA and AOB
archaea_ID_files = ''
bacteria_ID_files = ''

