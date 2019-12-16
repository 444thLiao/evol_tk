"""
test the correlation between phylogentic transition of habitats and the presence/absence of each ko.
"""
import pandas as pd
from os.path import *
import os
from tqdm import tqdm
from subprocess import check_call
import multiprocessing as mp
from glob import glob
from api_tools.itol_func import to_binary_shape

gene_presence_tab = "./protein_annotations/kegg_hmm_merged_info_e20.tab"
basic_habitat_txt = "./habitat_txt/habitat_general_cat.txt"
intree = "./bayestraits_habitat/over20p_bac120/over20p_bac120.formatted.newick"
odir = './bayestraits_habitat/over20p_bac120/gene_test_hmm'
habitat_mapping_dict = {"M":'1',
                        "N":'0',
                        "NM":'-'}
def run(cmd):
    check_call(cmd,shell=1,
               stdout=open("/dev/null",'w'))
    
habitat_text = open(basic_habitat_txt).read()
all_gids = [_.split('\t')[0] for _ in habitat_text.split('\n')]


gid2habitat = dict([_.split('\t') for _ in habitat_text.split('\n') if _])
genes_df = pd.read_csv(gene_presence_tab,sep='\t',index_col=0)


for ko in tqdm(genes_df.columns):
    gid2gene = genes_df.loc[:,ko].to_dict()
    gid2gene = {k:'1' if '_' in str(v) else '0' for k,v in gid2gene.items()}
    g_dir = join(odir,'each_gene',ko.split(':')[-1])
    if not exists(g_dir):
        os.makedirs(g_dir)
    metadata_txt = []
    for gid,v in gid2habitat.items():
        hv = habitat_mapping_dict[v]
        gene_v = str(int(gid2gene[gid]))
        if gene_v =='0':
            gene_v = '0'
        else:
            gene_v = '1'
        metadata_txt.append(f"{gid}\t{hv}\t{gene_v}")
    with open(join(g_dir,'metadata.txt'),'w') as f1:
        f1.write('\n'.join(metadata_txt))


# params1 = "3\n1\nRun" # Discrete: Dependant
# params2 = "2\n1\nRun" # Discrete: Independent
dependent_params = '/home-user/thliao/data/plancto/bayesTraits_genes_test/depend_params.txt'
independent_params = '/home-user/thliao/data/plancto/bayesTraits_genes_test/independ_params.txt'
exe_path = "/home-user/thliao/software/BayesTraitsV3.0.2-Linux/BayesTraitsV3"
cmds = []
for ko in tqdm(genes_df.columns):
    g_dir = join(odir,'each_gene',ko.split(':')[-1])
    mtext = join(g_dir,'metadata.txt')
    os.system(f"cp {mtext} {mtext.replace('.txt','_D.txt')}" )
    os.system(f"cp {mtext} {mtext.replace('.txt','_ID.txt')}" )
    cmds.append(f"{exe_path} {intree} {mtext.replace('.txt','_D.txt')} < {dependent_params}")
    cmds.append(f"{exe_path} {intree} {mtext.replace('.txt','_ID.txt')} < {independent_params}")



with mp.Pool(processes=50) as tp:
    r = list(tqdm(tp.imap(run,cmds),total=len(cmds)))

def read_log(infile,key='Tree No'):
    r_dict = {}
    rows = open(infile).readlines()
    idx,header = [(idx,header.strip('\n')) for idx,header in enumerate(rows) if header.startswith(key)][0]
    header = header.split('\t')
    for row in rows[idx+1:]:
        if row.strip('\n'):
            r_dict.update(dict(zip(header,row.strip('\n').split('\t'))))
    return r_dict
    
# collect results
ko2ID_rate = {}
ko2D_rate = {}
ko2lr = {}
for ko in tqdm(genes_df.columns):
    g_dir = join(odir,'each_gene',ko.split(':')[-1])
    Dresults = join(g_dir,'metadata_D.txt.Log.txt')
    IDresults = join(g_dir,'metadata_ID.txt.Log.txt')
    
    D_dict = read_log(Dresults)
    ID_dict = read_log(IDresults)
    ID_rate = {k:float(v) for k,v in ID_dict.items() if k.startswith('alpha') or k.startswith('beta')}
    D_rate = {k:float(v) for k,v in D_dict.items() if k.startswith('q')}
    ko2ID_rate[ko] = ID_rate
    ko2D_rate[ko] = D_rate
    l_D = float(D_dict['Lh'])
    l_ID = float(ID_dict['Lh'])
    lr = 2*(l_D-l_ID)
    ko2lr[ko] = lr
    
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests 
collect_ko2lr = {}
for ko,lr in ko2lr.items():
    if lr >=0:
        pval = 1 - chi2.cdf(lr,4)
        collect_ko2lr[ko] = pval
        
corrected_p = multipletests([_ for k,_ in collect_ko2lr.items()],
              method='bonferroni')
corrected_ko2p = dict(zip([k for k,_ in collect_ko2lr.items()],
                          corrected_p[1]))

significant_ko = [k for k,v in corrected_ko2p.items() if v<=0.05]
with open(join(odir,'significant_ko.list'),'w') as f1:
    f1.write('\n'.join(significant_ko))
kegg_info_tab = '/home-user/thliao/data/protein_db/kegg/ko_info.tab'
rows = [row.strip('\n') for row in open(kegg_info_tab) if row.split('\t')[0] in set(significant_ko)]
with open(join(odir,'significant_ko_info.tab'),'w') as f1:
    f1.write('\n'.join(rows))

subset_ko2D_rate = {ko:v for ko,v in ko2D_rate.items() if ko in significant_ko}
subset_df = pd.DataFrame.from_dict(subset_ko2D_rate,orient='index')
subset_df.to_csv(join(odir,"significant_ko_rate.tab"),sep='\t')

from sklearn.decomposition import PCA
pca = PCA()
pca_r = pca.fit_transform(subset_df.values)
pca_r = pd.DataFrame(pca_r,index=subset_df.index,columns=[f'PC{_+1}' for _ in range(pca_r.shape[1])])
pca_r.loc[:,'name'] = pca_r.index
g1 = pca_r.index[pca_r.loc[:,'PC1']>=20]
g2 = pca_r.index[pca_r.loc[:,'PC1']<20]

import plotly.express as px
fig = px.scatter(pca_r,x='PC1',y='PC2',text='name')
# fig.write_html('./test.html')


# fisher exact test
from scipy.stats import fisher_exact
ko2tab = {}
ko2odd_ratio = {}
distinct_ko = {}
for ko in tqdm(genes_df.columns):
    tab_num = []
    ko1_g = genes_df.index[genes_df.loc[:,ko].fillna('').str.contains('_')]
    ko0_g = genes_df.index[~genes_df.loc[:,ko].fillna('').str.contains('_')]
    
    tab_num.append([len([g 
                         for g in ko1_g 
                         if gid2habitat.get(g,'')==h]) 
                    for h in  ["M","N"]])
    tab_num.append([len([g for g in ko0_g if gid2habitat.get(g,'')==h]) for h in  ["M","N"]])
    ko2tab[ko]= tab_num
    oddratio,p = fisher_exact(tab_num,alternative="two-sided")
    distinct_ko[ko] = p
    ko2odd_ratio[ko] = oddratio
        
fe_corrected_p = multipletests([_ for k,_ in distinct_ko.items()],
              method='bonferroni')
fe_corrected_ko2p = dict(zip([k for k,_ in distinct_ko.items()],
                          fe_corrected_p[1]))  
set_fe_ko = set({k:v for k,v in fe_corrected_ko2p.items() if v<=0.05})
set_bst_ko = set(significant_ko)
intersect_ko = set_fe_ko.intersection(set_bst_ko)


# write out binary for visualization
t = sorted([(ko2odd_ratio[ko],ko) for ko in intersect_ko])[-30:]
t = [_[1] for _ in t]
tmp = genes_df.loc[:,t].to_dict(orient='index')
ssubset_g2ko = {gid:[] for gid in tmp}
{ssubset_g2ko[gid].append(ko) for gid,vdict in tmp.items() for ko,v in vdict.items() if '_' in str(v)}
text = to_binary_shape(ssubset_g2ko,info_name='test ko',omitted_other=True)
with open('./test_ko.binary.txt','w') as f1:
    f1.write(text)