
from glob import glob
import os,io
from collections import defaultdict
import pandas as pd
from os.path import *
from tqdm import tqdm
import pickle

Nif_list = ['K02591',
            'K02586',
            'K02588']
Nod_list = ['K14658',
            'K14659',
            'K14666',
            'K09695',
            'K09694',
            ]

Nif_list = ['K02584',
 'K02585',
 'K02586',
 'K02587',
 'K02588',
 'K02589',
 'K02590',
 'K02591',
 'K02592',
 'K02593',
 'K02594',
 'K02595',
 'K02596',
 'K02597',
 'K03737',
 'K03839',
 'K04488',
 'K15790']
Nod_list = [
 'K14657',
 'K14658',
 'K14659',
 'K14660',
 'K14661',
 'K14666',
 'K12546',
 'K18904',
 'K09694',
 'K09695',
 ]
# annotated_tab = "/home-user/jjtao/Rhizobiales/kegg_hmmsearch/e10/merged_result/merged_hmm_info.tab"
info_df = pd.read_csv("./nodnif_annotated_df.tsv",sep='\t',index_col=0)

gene2ko = {}
for _,row in info_df.iterrows():
    for genome,gene in row.items():
        if pd.isna(gene):
            continue
        for g in gene.split(','):
            gene2ko[g] = _
            
# sub_df = info_df.loc[Nif_list+Nod_list,:]
# sub_df.to_csv('./nodnif_annotated_df.tsv',sep='\t',index=1)


OG_df = pd.read_csv('./subset_og.tsv',sep='\t',index_col=1)
tmp_dir = './.tmp'
if exists(join(tmp_dir,'genome2gene_info')):
    genome2gene_info = pickle.load(open(join(tmp_dir,'genome2gene_info'), 'rb'))
    genome2order_tuple = pickle.load(open(join(tmp_dir,'genome2order_tuple'), 'rb'))
    genome2gene_info.pop('Oligotropha_carboxidovorans_OM5')
def get_genes(KO_list):
    subset_df = info_df.loc[KO_list,:]
    genes = [_ for v in subset_df.values for g in v for _ in str(g).split(',') if not pd.isna(g)]
    return genes

# get name of genes with annotated KO
Nod_genes = get_genes(Nod_list)
Nif_genes = get_genes(Nif_list)


## rename gbk (do only once)
import re
from ete3 import Tree
def refine_tree(tree_file,indir=None):
    text = open(tree_file).read()
    new_text = re.sub("\[[0-9]+\]",'',text)
    new_text = new_text.replace('OROOT','')
    t = Tree(new_text)
    right_ordered_names = t.get_leaf_names()
    if indir is not None and exists(indir):
        right_ordered_names = [rn
                               for rn in right_ordered_names
                               if exists(join(indir,f"{rn}.fna"))]
    return right_ordered_names
    
right_ordered_names = refine_tree("./genome_pruned.newick",
                                  './nif_neightbour100_genbank')

except_names = []
for rn in right_ordered_names:
    if not exists(join('./tmp_gbk', f'{rn}.gbk')):
        except_names.append(rn)
        
remap = {}  # from right name to wrong name
for _ in except_names:
    name = _.split('_')[-1]
    if len(glob(f'./tmp_gbk/*{name}.gbk'))!=1:
        print(_)
    else:
        remap[_] = glob(f'./tmp_gbk/*{name}.gbk')[0].split('/')[-1].replace('.gbk','')
         
remap['Bradyrhizobium_elkanii_TnphoA33'] = 'Bradyrhizobium_elkanii_TnphoA_33'
remap['Bradyrhizobium_sp_AS23_2'] = 'Bradyrhizobium_sp._AS23.2'
remap['Bradyrhizobium_sp_CNPSo_3424'] = 'Bradyrhizobium_sp._CNPSo_3424'
remap['Bradyrhizobium_sp_CNPSo_3426'] = 'Bradyrhizobium_sp._CNPSo_3426'
remap['Bradyrhizobium_sp_LMTR_3'] = 'Bradyrhizobium_sp._LMTR_3'
remap['Bradyrhizobium_sp_LVM_105'] = 'Bradyrhizobium_sp._LVM_105'
remap['Bradyrhizobium_sp_NAS80_1'] = 'Bradyrhizobium_sp._NAS80.1'
remap['Bradyrhizobium_sp_NAS96_2'] = 'Bradyrhizobium_sp._NAS96.2'
remap['Bradyrhizobium_sp_ORS_278_ORS278'] = 'Bradyrhizobium_sp._ORS_278'
remap['Bradyrhizobium_sp_ORS_375_ORS375'] = 'Bradyrhizobium_sp._ORS_375'
remap['Bradyrhizobium_sp_cf659_CF659'] = 'Bradyrhizobium_sp._cf659'
remap['Pseudolabrys_sp_GY_H'] = 'Pseudolabrys_sp._GY_H'
remap['Rhodopseudomonas_thermotolerans_JA576_NBRC_108863_KCTC_15144'] = 'Rhodopseudomonas_thermotolerans_JA576'
remap['Oligotropha_carboxidovorans_OM5'] = 'Oligotropha_carboxidovorans_OM5_ATCC_49405'
remap['Bradyrhizobium_sp_ORS_285'] = "Bradyrhizobium_sp._ORS_285"
remap['Rhodopseudomonas_thermotolerans_JA576_'] = 'Rhodopseudomonas_thermotolerans_JA576'


# assert
for rn in right_ordered_names:
    if not exists(join('./tmp_gbk', f'{rn}.gbk')):
        if rn not in remap:
            print(rn)

# fix names 
# for right,now_wrong_name in remap.items():
#     ori_name = join('./nif_neightbour100_genbank',f'{now_wrong_name}.fna')
#     new_name = join('./nif_neightbour100_genbank',f'{right}.fna')
#     if exists(ori_name):
#         os.system(f"mv {ori_name} {new_name}")

remap2 = {}
row = OG_df.iloc[0,:]
for col,v in row.items():
    v = v.split('|')[0]
    if v != col:
        remap2[v] = col
OG_df.columns = [_ if _ not in remap2 else remap2[_] for _ in OG_df.columns]

def count_number_contig(OG_df,genes,genome2gene_info,return_num=False):
    genome2num_contigs = defaultdict(set)
    for _ in genes:
        genome,gene = _.split('|')
        if genome not in genome2gene_info:
            genome = remap[genome]
        g2info = genome2gene_info[genome]
        contig = g2info[gene]['contig_name']
        genome2num_contigs[genome].add(contig)
    if return_num:
        genome2num_contigs = {k:len(v) for k,v in genome2num_contigs.items()}
    return genome2num_contigs


nod_genome2num_contigs = count_number_contig(OG_df,Nod_genes,genome2gene_info,return_num=True)
nif_genome2num_contigs = count_number_contig(OG_df,Nif_genes,genome2gene_info,return_num=True)

nod_genome2contigs = count_number_contig(OG_df,Nod_genes,genome2gene_info)
nif_genome2contigs = count_number_contig(OG_df,Nif_genes,genome2gene_info)

# nod and nif together
nod_nif_together_genomes = []
for genome,cset in nod_genome2contigs.items():
    cset2 = nif_genome2contigs.get(genome,set())
    if cset==cset2 and len(cset)==1:
        nod_nif_together_genomes.append(genome)
        
def process_locus_name_of_gbk2(row):
    infos = row.split(' ')
    infos = [_ for _ in infos if _]
    if len(infos)<=2:
        print(infos)
    try:
        length = infos[2]
    except:
        import pdb;pdb.set_trace()
    name = 'LOCUS' + ' '*7 + infos[1].split('_length')[0] + f' {length} bp  ' + '  DNA  linear  20-Jan-2020'
    return name+'\n'

def read_gbk(gbk_file):
    rows = open(gbk_file).readlines()
    rows = [process_locus_name_of_gbk2(_)
            if _.startswith('LOCUS') else _ 
            for _ in rows]
    rows = [process_locus_name_of_gbk2(_)
        if _.startswith('LOCUS') else _ 
        for _ in rows]
    records = SeqIO.parse(io.StringIO(''.join(rows)),format='genbank')
    records = list(records)
    return records

                
# extract_nifDKH
odir = './nif_neightbour100_genbank'
os.makedirs(odir,exist_ok=1)
neighbour = 50 * 1e3
for genome in tqdm([k for k,v in nif_genome2num_contigs.items() if v==1]):
    _nif_genes = [g for g in Nif_genes if remap.get(g.split('|')[0],g.split('|')[0]) == genome]
    f = f'./tmp_gbk/{genome}.gbk'
    if genome not in genome2gene_info:
        genome = remap[genome]
    g2info = genome2gene_info[genome]
    pos_list = []
    for gene in _nif_genes:
        gene = gene.split('|')[-1]
        pos1 = g2info[gene]['start']
        pos2 = g2info[gene]['end']
        contig = g2info[gene]['contig_name']
        pos_list+= [int(pos1),int(pos2)]
    if not pos_list:
        print(genome)
        continue
    if (max(pos_list) -min(pos_list)) >= 2*neighbour:
        continue
    
    start = int(min(pos_list) - neighbour )
    start = 0 if start <0 else start
    end = int(max(pos_list) + neighbour)
    ofile = join(odir,f'./{genome}.gbk')
    if exists(f) and not exists(ofile):
        genome_obj = [_ for _ in SeqIO.parse(f,format='genbank') if contig == _.id][0]
        new_genome = genome_obj[start:end]
        with open(ofile,'w') as f1:
            SeqIO.write([new_genome],f1,format='genbank')
    else:
        continue
for gbk in tqdm(glob(join(odir,'*.gbk'))):
    a = SeqIO.parse(gbk,format='genbank')
    with open(gbk.replace('.gbk','.fna'),'w') as f1:
            SeqIO.write(a,f1,format='fasta')
            
            
# stepwise blast
tree_file = "./nif_tree_pruned.newick"
tree_file = './brady_ortho354-345.txt'
text = open(tree_file).read()
new_text = re.sub("\[[0-9]+\]",'',text)
new_text = new_text.replace('OROOT','')
t = Tree(new_text)
right_names = t.get_leaf_names()

collect_rn = []
for rn in right_names:
    if exists(f"./nif_neightbour100_genbank/{rn}.fna"):
        collect_rn.append(f"./nif_neightbour100_genbank/{rn}.fna")

new_odir = join(odir,'genome_order_stepwise_align')
os.makedirs(new_odir,exist_ok=1)
for idx in tqdm(range(len(collect_rn))):
    next_idx = idx +1
    if next_idx == len(collect_rn)-1:
        break
    fna1 = collect_rn[idx]
    fna2 = collect_rn[idx+1]
    name1 = basename(fna1).replace('.fna','')
    name2 = basename(fna2).replace('.fna','')
    ofile = join(new_odir,name1+'_to_'+name2+'.blastout')
    if not exists(ofile):
        os.system(f"blastn -query {fna1} -subject {fna2} -evalue 1e-3 -outfmt 6 -out {ofile}")

# visulization
cmd = "/home-user/thliao/software/artemis/act "
for idx in list(range(len(collect_rn)))[:30]:
    next_idx = idx +1
    if next_idx == len(collect_rn)-1:
        break
    fna1 = collect_rn[idx]
    fna2 = collect_rn[idx+1]
    name1 = basename(fna1).replace('.fna','')
    name2 = basename(fna2).replace('.fna','')
    ofile = join(new_odir,name1+'_to_'+name2+'.blastout')
    cmd += f" {fna1} {ofile} "
cmd += f" {fna2} "
with open('./test.sh','w') as f1:
    f1.write(cmd)

# visulization via AliTV
text = ["""---
alignment:
  parameter:
    - --format=maf
    - --noytrim
    - --ambiguous=iupac
    - --gapped
    - --strand=both
  program: lastz
tree:   ./genome_pruned.newick"""]
text += ['genomes:']
for fna in collect_rn:
    name = basename(fna).replace('.fna','')
    text.append(f"  - name: {name}")
    text.append("    sequence_files:")
    text.append(f"      - {fna}")
with open('./genome_order_align.yml','w') as f1:
    f1.write('\n'.join(text))

# reformat the output fo AliTV
import json
json_obj = json.load(open('./genome_wide.json'))
json_obj['conf']['graphicalParameters']['canvasHeight'] = 20000
json_obj['conf']['graphicalParameters']['canvasWidth'] = 1000
json_obj['conf']['graphicalParameters']['fade'] = 0.8
json_obj['conf']['graphicalParameters']['karyoHeight'] = 30
json_obj['conf']['graphicalParameters']['fade'] = 0.8
## get name
name2seq_num = {}
for _seq,d in json_obj['data']['karyo']['chromosomes'].items():
    genome_name = d['genome_id']
    name2seq_num[genome_name] = _seq
## get order and stepwise relationship
tree_file = './nif_tree_pruned.newick'
tree_file = './brady_ortho354-345.txt'
text = open(tree_file).read()
new_text = re.sub("\[[0-9]+\]",'',text)
new_text = new_text.replace('OROOT','')
t = Tree(new_text)
right_names = t.get_leaf_names()
collect_rn = []
for rn in right_names:
    if exists(f"./nif_neightbour100_genbank/{rn}.fna") and rn in name2seq_num:
        collect_rn.append(rn)
json_obj['filters']['karyo']['genome_order'] = collect_rn
json_obj['filters']['karyo']['order'] = [name2seq_num[_]
                                         for _ in collect_rn]

# add features
json_obj["conf"]['features']['supportedFeatures'] = {"Nod":{
                    "color": "#ff006e",
                    "form": "rect",
                    "height": "30",
                    "visible": True
                },
                                             "Nif":{
                    "color": "#ffaf00",
                    "form": "rect",
                    "height": "30",
                    "visible": True
                }}


# generate features Table

nod_rows = []
nif_rows = []
for genome_name in tqdm(collect_rn):
    fna = f"./nif_neightbour100_genbank/{genome_name}.fna"
    genome = basename(fna).replace('.fna','')
    gbk = fna.replace('.fna','.gbk')
    record = SeqIO.read(gbk,format='genbank')
    c = record.id
    all_locus =[_.qualifiers['locus_tag'][0] for _ in record.features if _.type == 'CDS']
    _nod_genes = [g for g in Nod_genes 
                  if g.split('|')[-1] in all_locus]
    _nif_genes = [g for g in Nif_genes 
                  if g.split('|')[-1] in all_locus]
    if genome not in genome2gene_info:
        genome = remap[genome]
    g2info = genome2gene_info[genome]
    for _genes,_region in zip([_nod_genes,
                               _nif_genes],
                              [nod_rows,
                               nif_rows]):
        for fullgene in _genes:
            gene = fullgene.split('|')[-1]
            pos1 = g2info[gene]['start']
            pos2 = g2info[gene]['end']
            strand = g2info[gene]['strand']
            strand = '1' if strand == '+' else '-1'
            contig = g2info[gene]['contig_name']
            _region.append("\t".join(map(str,[genome_name,pos1,pos2,strand,gene2ko[fullgene]])))

nod_data = []            
for row in nod_rows:
    rows = row.split('\t')
    nod_data.append({"end":int(rows[2]),
                      "karyo":name2seq_num[rows[0]],
                      "start":int(rows[1]),
                      "name":"Nod"})
json_obj['data']['features']['Nod'] = nod_data
nif_data = []            
for row in nif_rows:
    rows = row.split('\t')
    nif_data.append({"end":int(rows[2]),
                      "karyo":name2seq_num[rows[0]],
                      "start":int(rows[1]),
                      "name":"Nif"})
json_obj['data']['features']['Nif'] = nif_data

obj2 = json.load(open('./genome_wide2.json'))
json_obj['data']['tree'] = obj2['data']['tree']
with open('./genome_wider_reset.json','w') as f1:
    json.dump(json_obj,f1)

# extract_fna
g2nod_region = defaultdict(list)
g2nif_region = defaultdict(list)
for genome in nod_nif_together_genomes:
    _nod_genes = [g for g in Nod_genes if remap.get(g.split('|')[0],g.split('|')[0]) == genome]
    _nif_genes = [g for g in Nif_genes if remap.get(g.split('|')[0],g.split('|')[0]) == genome]
    if genome not in genome2gene_info:
        genome = remap[genome]
    g2info = genome2gene_info[genome]
    for _genes,_region in zip([_nod_genes,
                               _nif_genes],
                              [g2nod_region,g2nif_region]):
        pos_list = []
        for gene in _genes:
            gene = gene.split('|')[-1]
            pos1 = g2info[gene]['start']
            pos2 = g2info[gene]['end']
            contig = g2info[gene]['contig_name']
            pos_list.append(pos1)
            pos_list.append(pos2)
            _region[genome].append(f'{contig}:{pos1}-{pos2}')
            

# too small to view
# gene distribution at the whole genome view
import plotly.graph_objs as go
from Bio import SeqIO
fig = go.Figure()
for g in g2nod_region:
    fig.add_scatter(x=[0,100],y=[g,g],mode='lines',line=dict(color='black',width=1),showlegend=False)
for _g2region in [g2nod_region,g2nif_region]:
    xs = []
    ys = []
    for genome,regions in tqdm(_g2region.items()):
        f = f'./tmp_gbk/{genome}.gbk'
        c = regions[0].split(':')[0]
        if exists(f):
            gb = [row for row in open(f).read().split('\n') if row.startswith('LOCUS') and c.split('.')[0] in row][0]
            idx = gb.split(' ').index('bp')
            num_length = int(gb.split(' ')[idx-1])
        else:
            continue
        for each_gene in regions:
            start,end = map(int,each_gene.split(':')[-1].split('-'))
            start,end = map(lambda x: x/num_length *100, [start,end])
            xs += [start,end,None]
            ys += [genome] *3
    fig.add_scatter(x=xs,y=ys,mode='lines+markers',line=dict(width=15))

fig.data[-2].marker.color = '#0000FF'
fig.data[-1].marker.color = '#FF0000'
fig.data[-2].name = 'Nod'
fig.data[-1].name = 'Nif'
fig.write_html('./test2.html')


# draw stack bar plot
def _tmp(x):
    x = x.split(':')[-1]
    s,e = map(int,x.split('-'))
    return abs(s-e)

fig = go.Figure()
for g in set(g2nod_region).union(set(g2nif_region)):
    regions1 = g2nod_region[g]
    regions2 = g2nif_region[g]
    _r = sorted(regions1+regions2,
                key=lambda x:int(x.split(':')[-1].split('-')[0]))
    total_length = sum([_tmp(_) for _ in _r])
    for gene in _r:
        length = _tmp(gene)
        if gene in regions1:
           fig.add_bar(x=[length/total_length *100],
                       y=[g],
                       orientation='h',
                       marker=dict(color='#0000FF')) 
        else:
           fig.add_bar(x=[length/total_length*100],
                       y=[g],
                       orientation='h',
                       marker=dict(color='#FF0000')) 
fig = fig.update_layout(barmode='stack')
fig.write_html('./test.html')