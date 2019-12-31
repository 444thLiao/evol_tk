import multiprocessing as mp
import os
from os.path import *
from subprocess import check_call

import pandas as pd
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
from glob import glob

def run(cmd):
    check_call(cmd, shell=1,
               stdout=open("/dev/null", 'w'))


def read_log(infile, key='Tree No'):
    r_dict = {}
    rows = open(infile).readlines()
    idx, header = [(idx, header.strip('\n')) for idx, header in enumerate(rows) if header.startswith(key)][0]
    header = header.split('\t')
    for row in rows[idx + 1:]:
        if row.strip('\n'):
            r_dict.update(dict(zip(header, row.strip('\n').split('\t'))))
    return r_dict


dependent_params = '/home-user/thliao/data/plancto/bayestraits_habitat/depend_params.txt'
independent_params = '/home-user/thliao/data/plancto/bayestraits_habitat/independ_params.txt'
exe_path = "/home-user/thliao/software/BayesTraitsV3.0.2-Linux/BayesTraitsV3"
gene_presence_tab = "./protein_annotations/hmmsearch_merged/merged_hmm_binary.tab"

indir = "./bayestraits_habitat/over20p_m2nm_v2/"

def batch_test_g(indir,odir=None):
    if odir is None:
        odir = join(indir,'gene_test')
    basic_habitat_txt = join(indir,"metadata.txt")
    intree = glob(join(indir,'*.newick'))[0]
    
    # read habitat metadata
    habitat_text = open(basic_habitat_txt).read()
    # read gene table 
    gid2habitat = dict([_.split('\t') for _ in habitat_text.split('\n')])
    genes_df = pd.read_csv(gene_presence_tab, sep='\t', index_col=0)

    # habitat mapping table...
    habitat_mapping_dict = {"M": '1',
                            "N": '0',
                            "NM": '-'}

    tqdm.write("Iterating genes to generating metadata.txt for each gene")
    for ko in tqdm(genes_df.columns):
        gid2gene = genes_df.loc[:, ko].to_dict()
        # below is for str....
        # gid2gene = {k: 1 if  is not None else 0 for k, _ in gid2gene.items()}
        g_dir = join(odir, 'each_gene', ko.split(':')[-1])
        if not exists(g_dir):
            os.makedirs(g_dir)
        metadata_txt = []
        for gid, v in gid2habitat.items():
            hv = habitat_mapping_dict[v]
            gene_v = str(int(gid2gene[gid]))
            metadata_txt.append(f"{gid}\t{hv}\t{gene_v}")
        with open(join(g_dir, 'metadata.txt'), 'w') as f1:
            f1.write('\n'.join(metadata_txt))

    tqdm.write("collecting metadata and generating file list")
    cmds = []
    for ko in tqdm(genes_df.columns):
        g_dir = join(odir, 'each_gene', ko.split(':')[-1])
        mtext = join(g_dir, 'metadata.txt')
        os.system(f"ln -sf ./{basename(mtext)} {mtext.replace('.txt', '_D')}")
        os.system(f"ln -sf ./{basename(mtext)} {mtext.replace('.txt', '_ID')}")
        if not exists(mtext.replace('.txt', '_D') + '.Log.txt'):
            cmds.append(f"{exe_path} {intree} {mtext.replace('.txt', '_D')} < {dependent_params}")
        if not exists(mtext.replace('.txt', '_ID') + '.Log.txt'):
            cmds.append(f"{exe_path} {intree} {mtext.replace('.txt', '_ID')} < {independent_params}")

    tqdm.write("run bayestraits!")
    with mp.Pool(processes=30) as tp:
        r = list(tqdm(tp.imap(run, cmds), total=len(cmds)))
    
from bin.transform.classify_kos import ko_classified_br,get_ko_infos,get_br_info
def parse_result(odir):
    # collect results
    ko2ID_rate = {}
    ko2D_rate = {}
    ko2lr = {}
    genes_df = pd.read_csv(gene_presence_tab, sep='\t', index_col=0)
    for ko in tqdm(genes_df.columns):
        g_dir = join(odir, 'each_gene', ko.split(':')[-1])
        Dresults = join(g_dir, 'metadata_D.Log.txt')
        IDresults = join(g_dir, 'metadata_ID.Log.txt')

        D_dict = read_log(Dresults)
        ID_dict = read_log(IDresults)
        ID_rate = {k: float(v) for k, v in ID_dict.items() if k.startswith('alpha') or k.startswith('beta')}
        D_rate = {k: float(v) for k, v in D_dict.items() if k.startswith('q')}
        ko2ID_rate[ko] = ID_rate
        ko2D_rate[ko] = D_rate
        l_D = float(D_dict['Lh'])
        l_ID = float(ID_dict['Lh'])
        lr = 2 * (l_D - l_ID)
        ko2lr[ko] = lr

    collect_ko2lr = {}
    for ko, lr in ko2lr.items():
        if lr >= 0:
            pval = 1-chi2.cdf(lr, 4)
            collect_ko2lr[ko] = pval
    list_ko = [k for k, _ in collect_ko2lr.items()]
    list_p = [collect_ko2lr[k] for k in list_ko]
    corrected_p = multipletests(list_p,
                                method='bonferroni')
    corrected_ko2p = dict(zip(list_ko,
                            corrected_p[1]))

    significant_ko = [k for k, v in corrected_ko2p.items() if v <= 0.05]
    with open(join(odir, 'significant_ko.list'), 'w') as f1:
        f1.write('\n'.join(significant_ko))

    br_kos = ko_classified_br(significant_ko)
    info_kos = get_ko_infos(significant_ko)
    infos = get_br_info(br_kos)
    df = pd.concat(infos,axis=0)
    df.loc[:,'des'] = [info_kos.get(_,'') for _ in df.index]
    
    subset_ko2D_rate = {ko: v for ko, v in ko2D_rate.items() if ko in significant_ko}
    subset_df = pd.DataFrame.from_dict(subset_ko2D_rate, orient='index')
    # gene 0 and N & gene 1 and M
    subset_df.loc[:,'sig rate M'] = subset_df.loc[:,['q21','q31','q34','q24']].sum(1) - subset_df.loc[:,['q12','q13','q43','q42']].sum(1)
    subset_df.to_csv(join(odir, "significant_ko_rate.tab"), sep='\t')
    basic_habitat_txt = join(indir,"metadata.txt")
    habitat_text = open(basic_habitat_txt).read()
    gid2habitat = dict([_.split('\t') for _ in habitat_text.split('\n')])
    M_g = set([g for g in gid2habitat if gid2habitat.get(g, '') == 'M'])
    N_g = set([g for g in gid2habitat if gid2habitat.get(g, '') == 'N'])

    all_over80 = list(subset_df.index[(subset_df.loc[:,['q24','q42']] >=80).all(1)])
    rate_sig_M = list(subset_df.index[(subset_df.loc[:,'sig rate M'] >=0)])
    
    set_ko = set(significant_ko)
    br2kos = ko_classified_br(set_ko)
    md_kos = ko_classified_module(set_ko)
    md2info = get_md_infos(md_kos.keys())
    ko2info = get_ko_infos(set_ko)
    collect_df = []

    for _subbr, k in tqdm(br2kos.items()):
        _subbr = {_subbr: k}
        df_list = get_br_info(_subbr)
        if not df_list:
            continue
        sub_df = pd.concat(df_list,axis=0,sort=False)
        sub_df.loc[:,'def'] = [ko2info[_] for _ in sub_df.index]
        collect_df.append(sub_df)

    tmp_df = pd.concat(collect_df, axis=0,sort=False)
    tmp_dict = {_:{'def':ko2info[_]}
     for _ in set_ko if _ not in tmp_df.index}
    tmp_df = tmp_df.append(pd.DataFrame.from_dict(tmp_dict,orient='index'),sort=False)
    for md,kos in md_kos.items():
        tmp_df.loc[kos,"module ID"] = md
        tmp_df.loc[kos,"module name"] = md2info[md]['name']
        tmp_df.loc[kos, "module all ko"] = md2info[md]['all ko']

    for ko,row in tmp_df.iterrows():
        ko1_g = genes_df.index[genes_df.loc[:, ko]==1]
        tmp_df.loc[ko,'Ratio (M) (%)'] = len(ko1_g.intersection(M_g))/len(M_g) *100
        tmp_df.loc[ko,'Ratio (N) (%)'] = len(ko1_g.intersection(N_g))/len(N_g) *100
        
    _all_over80 = [_ for _ in all_over80 if _ in tmp_df.index]
    _rate_sig_M = [_ for _ in rate_sig_M if _ in tmp_df.index]
    tmp_df.loc[_all_over80,'drive dual change'] = '+'
    tmp_df.loc[_rate_sig_M,'sig for M'] = subset_df.loc[_rate_sig_M,'sig rate M']
    tmp_df.loc[:,'LR'] = [ko2lr[_] for _ in tmp_df.index]
    tmp_df = tmp_df.sort_values(['top_br', 'module name',])
    tmp_df = tmp_df.reindex(columns=["LR",'Ratio (M) (%)','Ratio (N) (%)',"sig for M",'def',"module name",'module ID','module all ko',
                            'top','top_br','A','B','C','D',])
    tmp_df.to_csv(join(odir, f'significant_ko_info_only_bayes.tab'), sep='\t', index=1, index_label='K number')
    
    
    
# fisher exact test
from scipy.stats import fisher_exact

ko2tab = {}
ko2odd_ratio = {}
distinct_ko = {}
for ko in tqdm(genes_df.columns):
    tab_num = []
    ko1_g = genes_df.index[genes_df.loc[:, ko]==1]
    ko0_g = genes_df.index[genes_df.loc[:, ko]==0]

    tab_num.append([len([g
                         for g in ko1_g
                         if gid2habitat.get(g, '') == h])
                    for h in ["M", "N"]])
    tab_num.append([len([g for g in ko0_g if gid2habitat.get(g, '') == h]) for h in ["M", "N"]])
    ko2tab[ko] = tab_num
    oddratio, p = fisher_exact(tab_num, alternative="two-sided")
    distinct_ko[ko] = p
    ko2odd_ratio[ko] = oddratio

fe_corrected_p = multipletests([_ for k, _ in distinct_ko.items()],
                               method='bonferroni')
fe_corrected_ko2p = dict(zip([k for k, _ in distinct_ko.items()],
                             fe_corrected_p[1]))
set_fe_ko = set({k: v for k, v in fe_corrected_ko2p.items() if v <= 0.05})
set_bst_ko = set(significant_ko)
intersect_ko = set_fe_ko.intersection(set_bst_ko)

# classified and annotated with info
from bin.transform.classify_kos import get_br_info, ko_classified_br, get_ko_infos
import pandas as pd

M_g = set([g for g in gid2habitat if gid2habitat.get(g, '') == 'M'])
N_g = set([g for g in gid2habitat if gid2habitat.get(g, '') == 'N'])

all_over80 = list(subset_df.index[(subset_df.loc[:,['q24','q42']] >=80).all(1)])
rate_sig_M = list(subset_df.index[(subset_df.loc[:,'sig rate M'] >=0).all(1)])
for name,set_ko in [('only_bayes',set_bst_ko),('combine_fisher',intersect_ko)]:
    br2kos = ko_classified_br(set_ko)
    ko2info = get_ko_infos(set_ko)
    collect_df = []

    for _subbr, k in tqdm(br2kos.items()):
        _subbr = {_subbr: k}
        ko2brinfo = get_br_info(_subbr)
        if not ko2brinfo:
            continue
        ko2allinfo = {ko: {'def': info} 
                      for ko, info in ko2info.items() 
                      if ko in k}
        for ko, brinfo in ko2brinfo.items():
            ko2allinfo[ko].update(brinfo)
        sub_df = pd.DataFrame.from_dict(ko2allinfo, orient='index')
        collect_df.append(sub_df)

    tmp_df = pd.concat(collect_df, axis=0)
    tmp_dict = {_:{'def':ko2info[_]}
     for _ in set_ko if _ not in tmp_df.index}
    tmp_df = tmp_df.append(pd.DataFrame.from_dict(tmp_dict,orient='index'))
    tmp_df = tmp_df.sort_values(['top', 'A', 'B', 'C', 'D'])

    for ko,row in tmp_df.iterrows():
        ko1_g = genes_df.index[genes_df.loc[:, ko]==1]
        tmp_df.loc[ko,'M %'] = len(ko1_g.intersection(M_g))/len(M_g) *100
        tmp_df.loc[ko,'N %'] = len(ko1_g.intersection(N_g))/len(N_g) *100
        
        tmp_df.loc[ko,'M % over90Comple'] = len(ko1_g.intersection(over90_M_g))/len(over90_M_g) *100
        tmp_df.loc[ko,'N % over90Comple'] = len(ko1_g.intersection(over90_N_g))/len(over90_N_g) *100
        
        tmp_df.loc[ko,'Aquatic %'] = len(ko1_g.intersection(Aquatic_g))/len(Aquatic_g) *100
        tmp_df.loc[ko,'Non-Aquatic %'] = len(ko1_g.intersection(non_a_g))/len(non_a_g) *100
    _all_over80 = [_ for _ in all_over80 if _ in tmp_df.index]
    tmp_df.loc[_all_over80,'drive dual change'] = '+'
    tmp_df.to_csv(join(odir, f'significant_ko_info_{name}.tab'), sep='\t', index=1, index_label='K number')
    
# write out binary for visualization
t = sorted([(ko2odd_ratio[ko], ko) for ko in intersect_ko])[-30:]
t = [_[1] for _ in t]
tmp = genes_df.loc[:, t].to_dict(orient='index')
ssubset_g2ko = {gid: [] for gid in tmp}
{ssubset_g2ko[gid].append(ko) for gid, vdict in tmp.items() for ko, v in vdict.items() if '_' in str(v)}
text = to_binary_shape(ssubset_g2ko, info_name='test ko', omitted_other=True)
with open('./test_ko.binary.txt', 'w') as f1:
    f1.write(text)


# pca visualization
from sklearn.decomposition import PCA
import plotly.express as px
pca = PCA()
_t = subset_df.copy()
_t = _t.div(_t.sum(1),axis=0)
pca_r = pca.fit_transform(_t)
pca_r = pd.DataFrame(pca_r,
                     index=subset_df.index,
                     columns=[f'PC{_+1}' for _ in range(pca_r.shape[1])])
pca_r.loc[:,'name'] = pca_r.index
fig = px.scatter(pca_r,x='PC1',y='PC2',hover_name ='name')
fig.layout.xaxis.title.text = f'PC1 ({round(pca.explained_variance_ratio_[0]*100,2)} %)'
fig.layout.yaxis.title.text = f'PC2 ({round(pca.explained_variance_ratio_[1]*100,2)} %)'
fig.write_html(join(odir,'vis.html'))


# 
import itertools
from scipy.stats import pearsonr

pair2p = {}
for g1,g2 in tqdm(itertools.combinations(significant_ko,2)):
    g1_array = genes_df.loc[:,g1]
    g2_array = genes_df.loc[:,g2]
    r,p = pearsonr(g1_array,g2_array)
    pair2p[(g1,g2)] = p
    
list_pair = [k for k, _ in pair2p.items()]
list_p = [pair2p[k] for k in list_pair]
corrected_pair2p = multipletests(list_p,
                            method='bonferroni')
corrected_pair2p = dict(zip(list_pair,
                        corrected_pair2p[1]))

_subset_df = subset_df.copy()
_subset_df = _subset_df.div(_subset_df.sum(1),axis=0)
pair2p = {}
for g1,g2 in tqdm(itertools.combinations(significant_ko,2)):
    g1_array = _subset_df.loc[g1,:]
    g2_array = _subset_df.loc[g2,:]
    r,p = pearsonr(g1_array,g2_array)
    pair2p[(g1,g2)] = p
    
list_pair = [k for k, _ in pair2p.items()]
list_p = [pair2p[k] for k in list_pair]
corrected_pair2p = multipletests(list_p,
                            method='bonferroni')
corrected_pair2p = dict(zip(list_pair,
                        corrected_pair2p[1]))
rows = []
for (k1,k2),p in corrected_pair2p.items():
    if p<=0.05:
        rows.append(f'{k1}\t{k2}\t{p}')
with open('./bayestraits_habitat/phylo_pattern/network.tab','w') as f1:
    f1.write('\n'.join(rows))


from api_tools.itol_func import to_binary_shape
draw_data = pd.read_excel('./bayestraits_habitat/phylo_pattern/significant_ko_info_only_bayes.xlsx')
result_df = pd.read_csv('./protein_annotations/hmmsearch_merged/merged_hmm_binary.tab',sep='\t',index_col=0).T
all_gids = list(result_df.columns)

k2name = {}
for _,row in draw_data.iterrows():
    k2name[row['K number']] = row['def'].split(';')[0]

k2name['K20932'] = 'hzsA'
k2name['K20933'] = 'hzsB'
k2name['K20934'] = 'hzsC'
subset_df = draw_data.loc[~draw_data.label.isna()]
gb = subset_df.groupby('label')
for name,group_df in gb:
    kos = list(group_df.iloc[:,0])
    label = group_df.iloc[0,:]['label'].replace('/',' or ')
    shape = group_df.iloc[0,:]['shape']
    color = group_df.iloc[0,:]['color']

    id2ko = {}
    for gid in all_gids:
        remained_kos = set(result_df.index[result_df.loc[:,gid] > 0]).intersection(set(kos))
        id2ko[gid] = [k2name.get(_,_) for _ in list(remained_kos)]

    text = to_binary_shape(id2ko,
                           info_name=label,
                           manual_v=[k2name[_] for _ in kos],
                           info2style={k2name.get(k,k):{'shape':str(int(shape)),'color':color} for k in kos})
    with open(f'./bayestraits_habitat/phylo_pattern/itol_{label}.txt','w') as f1:
        f1.write(text)

