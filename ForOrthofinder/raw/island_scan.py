
from glob import glob
import os
from collections import defaultdict
import pandas as pd
from os.path import *
from tqdm import tqdm
import pickle

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
# sub_df = info_df.loc[Nif_list+Nod_list,:]
# sub_df.to_csv('./nodnif_annotated_df.tsv',sep='\t',index=1)


OG_df = pd.read_csv('./subset_og.tsv',sep='\t',index_col=1)
tmp_dir = './.tmp'
if exists(join(tmp_dir,'genome2gene_info')):
    genome2gene_info = pickle.load(open(join(tmp_dir,'genome2gene_info'), 'rb'))
    genome2order_tuple = pickle.load(open(join(tmp_dir,'genome2order_tuple'), 'rb'))

# test nod
subset_df = info_df.loc[Nod_list,:]
Nod_genes = [_ for v in subset_df.values for g in v for _ in str(g).split(',') if not pd.isna(g)]
# test nod
subset_df = info_df.loc[Nif_list,:]
Nif_genes = [_ for v in subset_df.values for g in v for _ in str(g).split(',') if not pd.isna(g)]

except_names = []
for _ in OG_df.columns:
    if not exists(join('./tmp_gbk', f'{_}.genbank')):
        print(_)
        except_names.append(_)

remap = {}
for _ in except_names:
    name = _.split('_')[-1]
    if len(glob(f'./tmp_gbk/*{name}.gbk'))!=1:
        # print(_)
        pass
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
remap['Bradyrhizobium_sp_ORS_285_ORS285'] = 'Bradyrhizobium_sp._ORS_285'
remap['Bradyrhizobium_sp_ORS_375_ORS375'] = 'Bradyrhizobium_sp._ORS_375'
remap['Bradyrhizobium_sp_cf659_CF659'] = 'Bradyrhizobium_sp._cf659'
remap['Pseudolabrys_sp_GY_H'] = 'Pseudolabrys_sp._GY_H'
remap['Rhodopseudomonas_thermotolerans_JA576_NBRC_108863_KCTC_15144'] = 'Rhodopseudomonas_thermotolerans_JA576'

def count_number_contig(OG_df,genes,genome2gene_info,return_num=False):
    genome2num_contigs = defaultdict(set)
    for _ in genes:
        genome,gene = _.split('|')
        if genome in remap:
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

# extract_fna
g2nod_region = {}
g2nif_region = {}
for genome in nod_nif_together_genomes:
    _nod_genes = [g for g in Nod_genes if remap.get(g.split('|')[0],g.split('|')[0]) == genome]
    _nif_genes = [g for g in Nif_genes if remap.get(g.split('|')[0],g.split('|')[0]) == genome]
    if genome in remap:
        genome = remap[genome]
    g2info = genome2gene_info[genome]
    for _genes,_region in zip([_nod_genes,
                   _nif_genes],[g2nod_region,g2nif_region]):
        pos_list = []
        for gene in _nod_genes+_nif_genes:
            gene = gene.split('|')[-1]
            pos1 = g2info[gene]['start']
            pos2 = g2info[gene]['end']
            contig = g2info[gene]['contig_name']
            pos_list.append(pos1)
            pos_list.append(pos2)
        if pos_list:
            _region[genome] = f'{contig}:{min(pos_list)}-{max(pos_list)}'