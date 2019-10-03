import pandas as pd
from collections import defaultdict
from sklearn.cluster import AgglomerativeClustering

tmp_path = '/home-user/thliao/data/metagenomes/update_0928_nitrification/confirmed_locus2info.tsv'
a = pd.read_csv(tmp_path,sep='\t',index_col=0)

            
genome2gene = defaultdict(list)
for _,row in a.iterrows():
    gid = row['locus_prefix']
    gene_name = row['Gene name(N metabolism)']
    genome2gene[gid].append(gene_name)
    
    

infile = 'test/query.tab'
annotate_tab = 'test/rep2duplicated.tab'
threshold = 0.2

def parse_df(infile,ofile,annotate_tab):
    total_df = pd.read_csv(infile,sep='\t',header=None)
    #remove_self_df = total_df.loc[~bool_index]
    mash_dis = [1 - float(mdist.split('/')[0])/float(mdist.split('/')[1]) 
                for mdist in total_df.iloc[:,-1].values]
    total_df.loc[:,'mash dist'] = mash_dis
    dist_df = total_df.pivot(index=0,columns=1,values='mash dist')
    dist_df.columns = [_.replace('./','') for _ in dist_df.columns]
    aggCluster = AgglomerativeClustering(affinity ='precomputed',distance_threshold =threshold,n_clusters=None,linkage ='single')
    cluster_labels = aggCluster.fit_predict(dist_df)
    
    group_dict = defaultdict(list)
    for clabel,name in zip(cluster_labels,dist_df.index):
        group_dict[clabel].append(name)
        
    group_dist_dict = defaultdict(list)
    for clabel,names in group_dict.items():
        group_dist_dict[clabel].append(pd.np.mean(dist_df.loc[names,names])[0])

    with open(annotate_tab,'w') as f1:
        f1.write('represent genome ID\tduplicated genome ID\tmash distance\n')
        for _,vlist in group_dict.items():
            if len(vlist) != 1:
                repID = vlist[0].split('/')[-1].replace('.fna','')
                for v in vlist[1:]:
                    altID = v.split('/')[-1].replace('.fna','')
                    dist = dist_df.loc[v,vlist[0]]
                    dist = str(round(dist,4))
                    f1.write(f'{repID}\t{altID}\t{dist}\n')
                    
    # get remained ID list
    remained_ID = [_[0].split('/')[-1].replace('.fna','')
                   for _ in group_dict.values()]
    dropped_ID = set([_ for vlist in group_dict.values() for _ in vlist]).difference(set(remained_ID))
    #return remained_ID
    
    
    
