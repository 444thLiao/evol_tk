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
import io
import multiprocessing as mp

odir = 'json_dump_v2'
indir = 'genome_protein_files_more'
og_tsv = f'./{indir}/OrthoFinder/Results_Oct01/Orthogroups/Orthogroups.tsv'
genome_info = './genome_info_full.xlsx'
MAG_annotate_file = '/home-user/thliao/data/metagenomes/update_0928_nitrification/confirmed_locus2info.tsv'
ref_file = '/home-user/thliao/project/nitrogen_cycle/nitrification/reference_genomes/outgroup and reference.xlsx'

# load all necessary data
with open(join(odir, 'ko2og.json'), 'r') as f1:
    ko2og = json.load(f1)
with open(join(odir, 'manually_blast_r.json'), 'r') as f1:
    manually_blast_r = json.load(f1)
with open(join(odir, 'ko2og2names.json'), 'r') as f1:
    ko2og2names = json.load(f1)

# load orthofinder csv, and rename it
og_df = pd.read_csv(og_tsv, sep='\t', low_memory=False, index_col=0)
def rename(x):
    if pd.isna(x):
        return x
    if ', ' not in x:
        return str(x).split(' ')[0]
    else:
        return ', '.join([str(_).split(' ')[0] for _ in x.split(', ')])
    
og_df = og_df.applymap(rename)
g_df = pd.read_excel(genome_info, index_col=0)
dropped_g_df = g_df.loc[g_df.loc[:,'used']=='no',:]
for gid in dropped_g_df.index:
    if gid in og_df.columns:
        og_df.drop(gid,axis=1,inplace=True)
        
# name2dirname
# for mag 
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

color_scheme = {'type':{'NOB': '#e41a1c', 'comammox': '#edc31d', 
                        'AOB': '#bad5b9', 'AOA': '#358f0f'},
                'phylum/class':{'Thaumarchaeota': '#358f0f',
                                'Nitrospirae': '#edc31d',
                                'Gammaproteobacteria': '#78fce0',
                                'Chloroflexi': '#e41a1c',
                                'Betaproteobacteria': '#956cb4',
                                'Alphaproteobacteria': '#8c613c',
                                'Actinobacteria':'#11FF11',
                                'Planctomycetes':'#FF66bb',
                                }

                }
phyla_color = {'Proteobacteria': '#358f0f',
                                'Firmicutes': '#edc31d',
                                'Actinobacteria': '#78fce0',
                                'Bacteroidetes': '#e41a1c'}


# rename
def rename_gene(each_og,ofile,extra={}):
    template_text = open(
        '/home-user/thliao/template_txt/labels_template.txt').read()
    id2org = generate_id2org(each_og, ofile)
    full_text = template_text[::]
    id2new_name = {}
    for id, org in id2org.items():
        if org in g_df.index:
            name = g_df.loc[org, 'genome name']
        else:
            name = id
        id2new_name[id] = name
    id2new_name.update(extra)
    return id2new_name


mag2info_file = '/home-user/thliao/project/nitrogen_cycle/nitrification/reference_genomes/mag2info.csv'
mag2info_df = pd.read_csv(mag2info_file,sep='\t',index_col=0)
mag2info_dict = mag2info_df.to_dict()


# data dependent transform
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
        
        
# ref or outgroup seq, additionally add to

def add_ref_seq(sub_df,used_records):
    used_ids = [str(_.id) for _ in used_records]
    id2info = {}
    record_need_dropped_ids = []
    for _,row in sub_df.iterrows():
        aa_id = row['AA accession'].strip()
        gene_name = row['gene name']
        seq = row['seq']
        info = row['phylum/class']
        if aa_id in used_ids:
            record_need_dropped_ids.append(aa_id)
            
        used_records.append(SeqIO.read(io.StringIO(f'>{aa_id}_{gene_name}\n{seq}'),format='fasta'))
        id2info[f'{aa_id}_{gene_name}'] = info
    final_records = [_ for _ in used_records if str(_.id) not in record_need_dropped_ids]    
    return final_records,id2info


def get_outgroup_info(sub_df):
    ID2name = {}
    ID2infos = {}
    for _,row in sub_df.iterrows():
        aa_id = row['AA accession'].strip()
        gene_name = row['gene name']
        org = row['org name']
        name = f'{aa_id}_{gene_name}'
        ID2name[name] = f'{org} ({gene_name})'
        if row['type'] == 'outgroup':
            ID2infos[name] = ['outgroup']
        else:
            ID2infos[name] = ['reference']
    info2style = {}
    info2style['outgroup'] = {'status':'0'}
    info2style['reference'] = {'status':'1'}
    return ID2infos,info2style,ID2name

def get_outgroup_info_phylum(sub_df,now_info2style):
    ID2infos = {}
    for _,row in sub_df.iterrows():
        aa_id = str(row['AA accession']).strip()
        gene_name = row['gene name']
        tax = str(row['phylum/class']).strip()
        name = f'{aa_id}_{gene_name}'
        ID2infos[name] = tax
    
    colors = px.colors.qualitative.Dark24 + px.colors.qualitative.Light24
    remained_colors = [c for c in colors if c not in now_info2style.values()]
    info2style = {}
    for v in set(ID2infos.values()):
        if v in now_info2style:
            pass
        elif v not in now_info2style and v in color_scheme:
            info2style.update({v:color_scheme[v]})
        elif v not in now_info2style and v not in color_scheme:
            one_color = remained_colors.pop(0)
            info2style.update({v:one_color})
    return ID2infos,info2style

def get_colors_general(ID2infos,now_info2style={}):
    colors = px.colors.qualitative.Dark24 + px.colors.qualitative.Light24
    remained_colors = [c for c in colors if c not in now_info2style.values()]
    info2style = {}
    for v in set(ID2infos.values()):
        if v in now_info2style:
            info2style[v] = now_info2style[v]
        elif v not in now_info2style and v in color_scheme:
            info2style.update({v:color_scheme[v]})
        elif v not in now_info2style and v not in color_scheme:
            one_color = remained_colors.pop(0)
            info2style.update({v:one_color})
    return ID2infos,info2style

# necessary for nxr and nar relative
# necessary for hao and hzo 
def refine_some_genes(fa_file,ko_name,no_dropped_ids=[]):
    removed_ids_files = glob(join('./manual_remove',ko_name+'*'))
    
    if removed_ids_files:
        removed_ids = [id
                       for f in removed_ids_files
                       for id in open(f).read().split('\n')]
    else:
        removed_ids = []
    
    records = [_ 
                for _ in SeqIO.parse(fa_file,format='fasta')
                if (_.id not in removed_ids) or (_.id in no_dropped_ids)]
    with open(fa_file.replace('.fa','.filterd.fa'),'w') as f1:
        SeqIO.write(records,f1,format='fasta-2line')
    print('refined ',fa_file)

outgroup_gene_names = {'K00370':['dms','tor'],
                       'K00371':['dms','tor'],
                       'K10535':['nrfA','_ONR'],
                       'K10944':['bmo'],
                       'K10945':['bmo'],
                       'K10946':['bmo']}

def process_ko(ko,og_list,final_odir,tree_exe='iqtree'):
    ref_df = pd.read_excel(ref_file,index_col=None)
    ref_df = ref_df.loc[ref_df.loc[:,'note']!='removed',:] 
    ofile = join(final_odir, ko+'.aln')
    ######
    sub_ref_df = ref_df.loc[ref_df.loc[:,'outgroup/ref for which KO']==ko,:]
    # for extra removed AOA
    sub_ref_df = sub_ref_df.loc[sub_ref_df.loc[:,'phylum/class']!='Thaumarchaeota',:]
    
    predir = dirname(dirname(og_tsv))

    # collect seq from Orthofinder and add outgroup sequence
    fa_files = [f'{predir}/Orthogroup_Sequences/{each_og}.fa' for each_og in og_list]
    used_records = [record for fa_file in fa_files for record in SeqIO.parse(fa_file,format='fasta')]
    final_records,ref_id2info = add_ref_seq(sub_ref_df,used_records)
    new_file = join(final_odir, ko+'.fa')
    with open(new_file,'w') as f1:
        SeqIO.write(final_records,f1,format='fasta-2line')
    
    # read manual remove directory, and remove the seqs want to remove        
    refine_some_genes(new_file,ko,
                      no_dropped_ids=list(ref_id2info.keys()))
    new_file = new_file.replace('.fa','.filterd.fa')
    final_ID_list = [_.id for _ in SeqIO.parse(new_file,format='fasta')]
    if not exists(ofile):
        check_call(f'mafft --maxiterate 1000 --genafpair --thread -1 {new_file} > {ofile}', shell=1)
        check_call(f"trimal -in {ofile} -out {ofile.replace('.aln','.trimal')} -automated1 -resoverlap 0.55 -seqoverlap 60", shell=1)
    if not exists(ofile.replace('.aln','.treefile')):
        #pass
        if tree_exe == 'iqtree':
            check_call(f"iqtree -nt AUTO -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre {final_odir}/{ko} -s {ofile.replace('.aln','.trimal')}",shell=1)
        else:
            n_file = ofile.replace('.aln','.treefile')
            check_call(f'FastTree {ofile} > {n_file}',shell=1)

    t = root_tree_with(ofile.replace('.aln','.treefile'),
                       gene_names=outgroup_gene_names.get(ko,[]),
                       format=0)
    renamed_tree(t,
                 outfile=ofile.replace('.aln','.newick'),
                 ascending=True)
    
    ## annotation part
    # outgroup and genome to habitat
    aa_ID2habitat = dict(zip( ['_'.join(map(lambda x:str(x).strip(),_)) 
                               for _ in ref_df.loc[:,['AA accession','gene name']].values],
                      ref_df.loc[:,'habitat (manual)'].str.strip()))
    org2habitat = dict(zip(g_df.index,
                           g_df.loc[:,'habitat (manual classify)']))
    
    id2org = generate_id2org(og_list, ofile)
    id2habitat = {id:org2habitat[org] 
                  for id,org in id2org.items()
                  if org in org2habitat}
    all_id2habitat = aa_ID2habitat.copy()
    all_id2habitat.update(id2habitat)
    # convert tree file output by iqtree and add internal node name
    
    # annotate the tree with bootstrap values as filled or unfilled dot.
    new_text = to_node_symbol(ofile.replace('.aln','.newick'))
    with open(join(final_odir,f'{ko}_node_bootstrap.txt'),'w') as f1:
        f1.write(new_text)
    
    ## annotate with type
    id2type = {k:g_df.loc[v,'type']
               for k,v in id2org.items()
               if v in g_df.index}
    id2info,info2col = get_colors_general(id2type,now_info2style= ref_id2info)
    write2colorstrip(id2info,final_odir,info2col,unique_id=ko,info_name='type')
    
    ## annotate with phylum/class as a color strip
    id2tax = {id: mag2info_dict['phylum/class'][gid]
                       for id,gid in id2org.items()
                       if gid in mag2info_dict['phylum/class']}
    _id2tax = {k:g_df.loc[v,'phylum/class']
               for k,v in id2org.items()
               if v in g_df.index}
    id2tax.update(_id2tax)
    id2tax = {k:v 
              for k,v in id2tax.items()
              if not pd.isna(v)}
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
    #id2tax = {k: v for k, v in id2tax.items() if v in Interested_tax}
    id2info,info2col = get_colors_general(id2tax,now_info2style= Interested_tax)
    _id2info,_info2col = get_outgroup_info_phylum(sub_ref_df,info2col)
    id2info.update(_id2info)
    info2col.update(_info2col)
    id2info = {k: v for k, v in id2info.items() if v in Interested_tax}
    write2colorstrip(id2info,final_odir,Interested_tax,unique_id=ko,info_name='phylum/class')
    
    # annotate branch color as the same as phylum/class with tree
    
    write2colorbranch_clade(id2info,
                            final_odir,
                            Interested_tax,
                            treefile=ofile.replace('.aln','.newick'),
                            unique_id=ko,
                            info_name='branch_color',
                            no_legend=False)
    # get new add sequence infomation and annotate it
    ID2infos,_,extra_id2name = get_outgroup_info(sub_ref_df)
    locus2new_name = rename_gene(og_list,ofile,extra=extra_id2name)
    full_text = to_label(locus2new_name)
    with open(join(final_odir, f'{ko}_label.txt'), 'w') as f1:
        f1.write(full_text)
        
    ID2add_type = {ID:ID2infos[ID][0] 
                   if ID2infos.get(ID,[]) else 'MAGs' 
                   for ID in final_ID_list}
    ID2add_type.update({ID:'reference genome' 
                   for ID in final_ID_list
                   if id2org.get(ID,'nan') in g_df.index})
    template_text = to_color_strip(
    ID2add_type, {"reference":'#FF0000',
                                                        "reference genome":'#FF0000',
                                                        "outgroup":'#00FF00',
                                                        "MAGs":'#0000FF'}, info_name='gene name')
    with open(join(final_odir, f'{ko}_from_strip.txt'), 'w') as f1:
        f1.write(template_text)
        
    ## annotate with habitat
    new_locus2habit = {id: mag2info_dict['habitat (manual)'][gid]
                       for id,gid in id2org.items()
                       if gid in mag2info_dict['habitat (manual)']}
    id2info,info2col = get_colors_general(new_locus2habit)
    all_id2habitat.update(id2info)
    all_id2habitat = {k:str(v).strip() for k,v in all_id2habitat.items()
                      if k in final_ID_list}
    all_id2habitat = {k:'terrestrial' if v =='soil' else v
                      for k,v in all_id2habitat.items()}
    all_id2habitat = {k:'freshwater' if v =='fresh water' else v
                      for k,v in all_id2habitat.items()}
    
    habitats = list(set(all_id2habitat.values()))
    habitat_colros = {'marine': '#0011FF',
                  'marine(symbiont)': '#0011FF',
                  'terrestrial': '#f99806',
                  'terrestrial(symbiont)': '#f99806',
                  'freshwater': '#88bbFF',
                  'waste water': '#50AF6E',
                  'artifical system':'#FA3705'
                  }
    matrix_text = to_matrix_shape(all_id2habitat,'habitat',color=habitat_colros)
    
    with open(join(final_odir,f'{ko}_habitat.txt'),'w') as f1:
        f1.write(matrix_text)
    ## annotate with gene name 
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
    template_text = to_color_strip(
        id2info, info2col, info_name='gene name')
    with open(join(final_odir, f'{ko}_gene_name_strip.txt'), 'w') as f1:
        f1.write(template_text)
    #habitat_filewrite2colorstrip(id2info,final_odir,info2col,unique_id=ko,info_name='habitat')
    
    
final_odir = join('./align_v3', 'complete_ko')
os.makedirs(final_odir, exist_ok=1)
params_list = []
for ko,og_list in ko2og.items():
    og_list = ko2og[ko]
    og_list = og_list[::]
    if ko == 'K10944':
        og_list.remove('OG0006386')  # arachaea amoA
        #og_list.append('OG0001887')  # amoA of Heterotrophic nitrification
        
    params_list.append((ko,og_list,final_odir))
    
def run_c(x):
    process_ko(*x)

with mp.Pool(processes=10) as tp:
    tp.map(run_c,params_list)



final_odir = join('./align_v3', 'trees_removed_AOA')
os.makedirs(final_odir, exist_ok=1)
params_list = []
for ko,og_list in ko2og.items():
    if ko.startswith('K109'):
        og_list = ko2og[ko]
        og_list = og_list[::]
        if ko == 'K10944':
            og_list.remove('OG0006386')  # arachaea amoA
            #og_list.append('OG0001887')  # amoA of Heterotrophic nitrification
            
        params_list.append((ko,og_list,final_odir))
    
def run_c(x):
    process_ko(*x)

with mp.Pool(processes=10) as tp:
    tp.map(run_c,params_list)
    

# get_all_genome_name
org_names = []
final_odir = join('./align_v3', 'complete_ko')
os.makedirs(final_odir, exist_ok=1)
params_list = []
for ko,og_list in ko2og.items():
    og_list = ko2og[ko]
    og_list = og_list[::]
    if ko == 'K10944':
        og_list.remove('OG0006386')  # arachaea amoA
        #og_list.append('OG0001887')  # amoA of Heterotrophic nitrification
        
    ofile = join(final_odir, ko+'.aln')
    id2org = generate_id2org(og_list, ofile)
    org_names += list(set(id2org.values()))
fa_files = []
for _ in set(org_names):
    flist = glob('./genome_protein_files_more/'+_+'*')
    if not flist:
        print(_)
    else:
        fa_files.append(flist[0])
        
from tqdm import tqdm
cdd_odir = join('./genome_protein_files_more','cogRbp_anno')
os.makedirs(cdd_odir,exist_ok=1)
for fa in tqdm(fa_files):
    anno = join(cdd_odir,basename(fa).replace('.faa','.cogrbp'))
    cmd_template = f"blastp -query {fa} -db /home-user/sswang/resource/db/CogRbp/CogRbp -outfmt 6 -num_threads 50 -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 > {anno}" 
    check_call(cmd_template,shell=1)
_sub_df =g_df.loc[:,['phylum/class','habitat (manual classify)']]
_sub_df.columns = ['phylum/class', 'habitat (manual)']
new_df = pd.concat([mag2info_df,_sub_df])
new_df.loc[:,'habitat (manual)'] = [_.strip() for _ in new_df.loc[:,'habitat (manual)']]
new_df.index = [_.split('.')[0] for _ in new_df.index]
new_df.to_csv('./genome_protein_files_more/extract_cog_aln/info.csv',sep='\t',index=1)


odir = './genome_protein_files_more/extract_cog_aln/'
all_ids = [_.id for _ in SeqIO.parse(join(odir,'concat_aln.aln'),format='fasta')]

id2tax = {id:new_df.loc[id,'phylum/class'] if id in new_df.index else new_df.loc[name2dirname[id].split('.')[0],'phylum/class'] for id in all_ids}
id2tax = {k:v 
            for k,v in id2tax.items()
            if not pd.isna(v)}
id2tax = {k:'CPR' if 'candidatus' in str(v).lower() or 'candidate' in str(v) else v
            for k,v in id2tax.items()}
id2info,info2col = get_colors_general(id2tax,now_info2style= color_scheme['phylum/class'])
write2colorstrip(id2info,odir,info2col,unique_id='concat_no_par',info_name='phylum/class')
root_with ='GCA_002345125'
    
write2colorbranch_clade(id2info,
                            odir,
                            info2col,
                            treefile=join(odir,'concat_no_par.newick'),
                            unique_id='concat_no_par',
                            info_name='branch_color',
                            no_legend=True)    

new_text = to_node_symbol(join(odir,'concat_no_par.newick'))
with open(join(odir,'concat_no_par_node_bootstrap.txt'),'w') as f1:
    f1.write(new_text)

    