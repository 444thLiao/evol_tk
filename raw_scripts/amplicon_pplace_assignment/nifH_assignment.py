
import copy
import os
import re
from collections import defaultdict
from glob import glob
from os.path import *
from subprocess import check_call
from ete3 import Tree
from Bio import SeqIO
from Bio.Data import IUPACData
from Bio.Seq import Seq
from tqdm import tqdm
import pandas as pd
import json
from bin.format_newick import renamed_tree
import numpy as np
tmp_d = IUPACData.ambiguous_dna_values

# r2seq = {_.id:_ for _ in SeqIO.parse('/mnt/maple/sswang/Brady_zong/Brady/results/gene_tree/official-3/pep/nifH.fas','fasta')}
# r2gene = {}
# for r in tqdm(r2seq):
#     g = [f"/home-user/sswang/project/Brady/data/Bradyrhizobiaceae/cds/{r}.gene",
#          f"/mnt/maple/sswang/Brady_zong/Brady/data/zong/cds/{r}.gene"]
#     g = [_ for _ in g if exists(_)]
#     if not g:
#        # continue
#         print(r)
#     records = list(SeqIO.parse(g[0],format='fasta'))
#     for gene in records:
#         if str(gene.seq.translate())[:-1]==str(r2seq[r].seq):
#             r2gene[r] = gene



# with open('/home-user/lling/5_Amplicon/nifH_tree/ref_nifH.fasta','w') as f1:
#     SeqIO.write(r2gene.values(),f1,'fasta-2line')
# with open('/home-user/thliao/tmp_nifH','w') as f1:
#     for r in tqdm(r2seq):
#         g = [f"/home-user/sswang/project/Brady/data/Bradyrhizobiaceae/cds/{r}.gene",
#             f"/mnt/maple/sswang/Brady_zong/Brady/data/zong/cds/{r}.gene"]
#         g = [_ for _ in g if exists(_)]
#         if not g:
#         # continue
#             f1.write(r+'\n')

os.chdir('/home-user/lling/5_Amplicon/PolF_Soil_seq_data/')

ref_newick = "/home-user/lling/5_Amplicon/built_nifH_db/merged_nifH.newick"
tre = Tree(ref_newick)
tre.resolve_polytomy()
tre = renamed_tree(tre)
# with open('./ref_nifH_resolved.newick','w') as f1:
tre.write(outfile='./ref_nifH_resolved.newick',format=3)
ref_newick = '/home-user/lling/5_Amplicon/PolF_Soil_seq_data/ref_nifH_resolved.newick'

n2plength = {n.name:n.dist for n in tre.traverse()}
all_dist = list(n2plength.values())
np.mean(all_dist),np.max(all_dist)
percentile95 = np.percentile(all_dist,95)

## trim the alignment
ref_aln = "/home-user/lling/5_Amplicon/built_nifH_db/merged_nifH.aln"
pre_sequences = list(SeqIO.parse(ref_aln,'fasta'))

primer_set = {"f": "TGCGAYCCSAARGCBGACTC",
              "r": "ATSGCCATCATYTCRCCGGA"}

def get_region(seq, f, r,include_primer=False):
    f_len, r_len = len(f), len(r)
    f = ''.join([l if len(tmp_d.get(l, l)) ==
                 1 else f"[{tmp_d.get(l,l)}]" for l in f])
    r = str(Seq(r).reverse_complement())
    f = ''.join([l if len(tmp_d.get(l, l)) ==
                 1 else f"[{tmp_d.get(l,l)}]" for l in f])
    r = ''.join([l if len(tmp_d.get(l, l)) ==
                 1 else f"[{tmp_d.get(l,l)}]" for l in r])
    if include_primer:
        matched_seq = re.search(f'({f}[ATCG]*{r})', seq)
    else:
        matched_seq = re.search(f'{f}([ATCG]*){r}', seq)
    if matched_seq is None:
        return None
    # return the sequence removing the primer
    return matched_seq.groups()[0]

def trim_aln_with_seq(ref_seq, regional_seq):
    # 0-coordinated, thus left closed and right closed
    aft2ori_idx = {}
    ori_i = aft_i = 0
    for s in ref_seq:
        if s.upper() not in 'ACGTN':
            ori_i += 1
        else:
            aft_i += 1
            ori_i += 1
        aft2ori_idx[aft_i] = ori_i
    ori_idx = list(range(0, len(ref_seq)))
    nogap_refseq = ''.join([_.upper()
                            for _ in ref_seq if _.upper() in 'ACGTN'])
    aft_idx = list(range(0, len(nogap_refseq)))

    idx = nogap_refseq.find(regional_seq.upper())
    if idx == -1:
        raise IOError
    start = aft2ori_idx[idx]
    end = aft2ori_idx[idx+len(regional_seq)]
    return start, end

matched_ = defaultdict(int)
for r in tqdm(pre_sequences):
    pos = [_ for _,bp in enumerate(r.seq) if bp !='-']
    new_s = [bp for _,bp in enumerate(r.seq) if bp !='-']
    pos_d = dict(zip(range(len(new_s)),pos))
    region = get_region(''.join(new_s),
                        primer_set['f'],
                        primer_set['r'],
                        include_primer=True)
    if not region:
        continue    
    s, e = trim_aln_with_seq(''.join(new_s), region)
    ns = pos_d[s]
    ne = pos_d[e]
    matched_[(ns,ne)] +=1

### manual check whether it has only one s,e
s,e = list(matched_.keys())[0]
print(s,e)
trimmed_seqs = []
for r in pre_sequences:
    sub_r = copy.deepcopy(r)[s:e]
    sub_r.seq = Seq(''.join([_ for _ in sub_r.seq if _!='-']))
    trimmed_seqs.append(sub_r)
with open('./ref_trimmed_nifH.fasta','w') as f1:
    SeqIO.write(trimmed_seqs,f1,'fasta-2line')
os.system(f"mafft --auto ./ref_trimmed_nifH.fasta > ./ref_trimmed_nifH.aln")
ref_phy = '/home-user/lling/5_Amplicon/PolF_Soil_seq_data/final_reference.phy'
records = list(SeqIO.parse('ref_trimmed_nifH.aln','fasta'))
n1,n2 = len(records),len(records[0].seq)
with open(ref_phy,'w') as f1:
    f1.write(f"{n1} {n2}\n")
    for r in records:
        f1.write(f"{r.id}{' '*10}{r.seq}\n")

## validations
os.system(f"cat ./ref_trimmed_nifH.fasta ../PolF_Soil_seq_data/Brady_nifH.fasta > ./ALL_together.fasta; mafft --auto ./ALL_together.fasta > ./ALL_together.aln")
os.system(f"FastTreeMP ./ALL_together.aln > ./ALL_together.newick")
os.system(f"mkdir trees; iqtree -nt 15 -m MFP -redo -mset GTR,HKY85,K80 -s ./ALL_together.aln -pre ./trees/ALL_together -fast ")
########
os.chdir('/home-user/lling/5_Amplicon/PolF_Soil_seq_data/assignments')
infa = '/home-user/lling/5_Amplicon/PolF_Soil_seq_data/output_dir/dada2_output/rep.fa'
if exists('./papara_log.aln'):
    os.system(f"rm ./papara_*.aln")
cmd1 = f"papara -t {ref_newick} -s {ref_phy} -q {infa} -r -n aln -j 30"
check_call(cmd1,shell=1)
cmd2 = f"epa-ng --split {ref_phy} papara_alignment.aln --redo"
check_call(cmd2,shell=1)
cmd3 = f"epa-ng --ref-msa ./reference.fasta --tree {ref_newick} --query ./query.fasta -T 30 -w ./  --model 'GTR+F+R10' --redo"
check_call(cmd3,shell=1)
# ./reference.fasta and query.fasta are generated


from for_software.for_EPA.parse_jplace import parse_tree_with_edges

guppy_exe = "/home-user/thliao/download/pplacer-Linux-v1.1.alpha19/guppy "
jplace = './epa_result.jplace'
cmd = f"{guppy_exe} to_csv {jplace} > {dirname(jplace)}/jplace.csv"
check_call(cmd,shell=1)
cmd = f"guppy edpl {jplace} --csv -o {dirname(jplace)}/edpl.csv"
check_call(cmd,shell=1)
edpl_df = pd.read_csv(f"{dirname(jplace)}/edpl.csv",index_col=0,header=None)
df = pd.read_csv(f"{dirname(jplace)}/jplace.csv")
df = df.sort_values('like_weight_ratio',ascending=False).groupby('name').head(1)
df.index = df['name']
df = df.reindex(edpl_df.index)
df.loc[:,'EDPL'] =  edpl_df[1]
df = df.reindex(columns=['edge_num', 'like_weight_ratio','distal_length', 'pendant_length','EDPL'])
obj = json.load(open(jplace))
used_tree = obj["tree"]
tree,node_name2edge_num = parse_tree_with_edges(used_tree)
e2n = {v:k for k,v in node_name2edge_num.items()}
df.loc[:,'node'] = [e2n[_] for _ in df['edge_num']]
df.loc[:,'parental length'] =  [n2plength[_] for _ in df['node']]
final_df = df.loc[(df['pendant_length']<=percentile95) & (df['EDPL']<=0.1),:] 
validated_seqs = list(final_df.index)

### assignments
# def classification_criteria(tree):
#     name2node = {_.name:_ for _ in tree.traverse()}
#     all_names = list(name2node)
#     pb_group = list(set(all_names).difference([_.name for _ in name2node['I132_S0'].traverse()]))
#     kakadu_group = ['Bradyrhizobium_sp_ARR65','Bradyrhizobium_axai_SZCCHNG1001']
    
#     jicamae_group = ['Bradyrhizobium_jicamae_PAC68']
#     for target_n in ['I250_S0','I169_S0']:
#         jicamae_group+=list([_.name for _ in name2node[target_n].traverse()])
    
#     japonicum_group = []
#     for target_n in ['I369_S0','I248_S0','I194_S1','I164_S0','I414_S0',]:
#         japonicum_group+=list([_.name for _ in name2node[target_n].traverse()])
#     elkanii_group = [_.name for _ in name2node['I165_S1'].traverse()]
#     elkanii_group = list(set(elkanii_group).difference(jicamae_group+japonicum_group))
#     criteria = {"pb":pb_group,
#                 "jicamae":jicamae_group,
#                 "japonicum":japonicum_group,
#                 'elkanii':elkanii_group}
#     return criteria
def classification_criteria(tree):
    name2node = {_.name:_ for _ in tree.traverse()}
    all_names = list(name2node)
    
    FL1_group = []
    for target_n in ['I164_S0']:
        FL1_group+=list([_.name for _ in name2node[target_n].traverse()])
    _a = list([_.name for _ in name2node['I211_S0'].traverse()])
    FL1_group = set(FL1_group).difference(set(_a))
    
    pb_group = list(set(all_names).difference([_.name for _ in name2node['I132_S0'].traverse()]))
    
    FL2_group = []
    for target_n in ['I161_S0']:
        FL2_group+=list([_.name for _ in name2node[target_n].traverse()])
    
    Sym_group =  list([_.name for _ in name2node['I154_S0'].traverse()])
    Sym_group = set(Sym_group).difference(set(FL1_group))

    criteria = {"pb":pb_group,
                "FL1":FL1_group,
                "FL2":FL2_group,
                'Sym':Sym_group}
    return criteria

criteria = classification_criteria(tre)
seq2type = {}
for type_name, nodes in criteria.items():
    sub_df = final_df.loc[final_df['node'].isin(nodes)]
    names = list(sub_df.index)
    seq2type.update({n: type_name for n in names})
# from api_tools.itol_func import *
# text = to_color_strip({node_name:n for n,ns in criteria.items() for node_name in ns},
#                       {"pb":'#9dd0f0',
#                        "FL2":'#ffd29a',
#                        "FL1":'#b4a7d6',
#                        'Sym':'#9addba'})
# with open('./tmp.txt','w') as f1:f1.write(text)
# text = to_color_strip({node_name:n for n,ns in criteria.items() for node_name in ns},
#                       {"pb":'#9dd0f0',
#                        "jicamae":'#ffd29a',
#                        "japonicum":'#b4a7d6',
#                        'elkanii':'#9addba'})
# with open('./tmp.txt','w') as f1:f1.write(text)


count_tab = f"/home-user/lling/5_Amplicon/PolF_Soil_seq_data/output_dir/dada2_output/profiling.csv"
count_df = pd.read_csv(count_tab, sep='\t', index_col=0)
count_df = count_df.T
count_df = count_df.loc[validated_seqs,:]
sample2total_reads = count_df.sum(0)
ratio_df = count_df/sample2total_reads * 100
sub_ratio_df = ratio_df.loc[ratio_df.index.isin(seq2type), :]
sub_ratio_df.loc[:, 'type'] = [seq2type[s] for s in sub_ratio_df.index]
type2ratio = sub_ratio_df.groupby('type').sum().T
type2ratio.columns = [f"{_} (%)" for _ in type2ratio.columns]
type2ratio.loc[:, 'reads inferred as nifH'] = count_df.sum(0).reindex(type2ratio.index)


sample2stats = pd.read_csv(dirname(count_tab)+'/profiling_stats.csv',sep='\t',index_col=0)
sample2stats = sample2stats.reindex(type2ratio.index)
sample2stats = sample2stats.reindex(columns=['input','filtered','denoised','merged','non-chimeric'])
sample2stats.columns=['input','After QC','dada2-denoised', 'dada2-merged','dada2-ASV']
type2ratio = type2ratio.join(sample2stats)
type2ratio.to_csv('./nifH_abundances.tsv',sep='\t',index=1)


#### TODO
# implemtn a taxonomy table and use gappa/guppy to assign the query sequences into designed groups

template = open('/mnt/home-backup/thliao/cyano/analysis/nfixer/dataset/refS1/gene2tax.txt').read().strip().split('\n')
refS1 = [_.id for _ in SeqIO.parse('/home-user/lling/5_Amplicon/built_nifH_db/refS1_nifH_nuc.fasta','fasta')]
missing_taxa = []
for r in tre.get_leaf_names():
    if r not in refS1:
        missing_taxa.append(r)
        
name2group = {}    
for group,n_list in criteria.items():
    for n in n_list:
        if n in missing_taxa:
            name2group[n] = group
        
from bin.ncbi_convertor.toolkit import tax2tax_info
from ete3 import NCBITaxa
ncbi = NCBITaxa()
name2tid = {}
for rawname in tqdm(missing_taxa):
    name = ' '.join(rawname.split('_',2)[:2])
    tid = ncbi.get_name_translator([name])
    if not tid:
        name = ' '.join(rawname.split('_',2)[:1])
        tid = ncbi.get_name_translator([name])
        name2tid[rawname] = tid[name][0]
    else:
        name2tid[rawname] = tid[name][0]

from bin.ncbi_convertor.toolkit import tax2tax_info
all_tids = set(name2tid.values())
i = ['superkingdom','phylum','order','family','genus']
for name,tid in name2tid.items():
    if str(tid)=='nan':continue
    info = tax2tax_info(tid)
    names = [info.get(i,'') for i in info]
    if name in name2group:
        template.append(f"{name}\t{';'.join(names+[name2group[name]])}")
with open('./ref2tax.txt','w') as f1:
    f1.write('\n'.join(template))
    
cmd = "gappa examine assign --jplace-path ./epa_result.jplace --taxon-file ./ref2tax.txt --allow-file-overwriting  --per-query-results  --ranks-string 'superkingdom|phylum|class|order|family|genus|species|group'  --distant-label  --sub-taxopath 'Bacteria;Proteobacteria'  --out-dir ./gappa "    
os.system(cmd)

ga_df = pd.read_csv('per_query.tsv',sep='\t')
ga_df = ga_df.sort_values('afract',ascending=False).groupby('name').head(1)
ga_df.index = ga_df['name']
ga_df = ga_df.drop(columns='name')
sub_ga_df = ga_df.loc[ga_df['taxopath'].str.contains('Bacteria;Proteobacteria;Alphaproteobacteria'),:]
    