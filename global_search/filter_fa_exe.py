from Bio import SeqIO,Entrez
from os.path import *
from tqdm import tqdm
import os
import numpy as np
from subprocess import check_call
import sys
from Bio import Entrez
from bin.ncbi_convert import edl
import io
import pandas as pd
def run(cmd):
    check_call(cmd,shell=True,stdout=open('/dev/null'))
kofam_scan = '/home-user/thliao/software/kofamscan/exec_annotation'


gene_info = {'kegg': {'nxrA': 'K00370',
                      'nxrB': 'K00371',
                      'hao': 'K10535',
                      'amoA': 'K10944',
                      'amoB': 'K10945',
                      'amoC': 'K10946'},
             'TIGFAM': {'nxrA': '',
                        'nxrB': '',
                        'hao': 'TIGR01703',
                        'amoA': 'TIGR03080',
                        'amoB': 'TIGR03079',
                        'amoC': 'TIGR03078'}}

def filter_by_relative_pos(full_df):
    pid2nuccore_ids = dict(zip(full_df.index,full_df.loc[:,'nuccore ID']))
    nuccore_ids = [v for k,v in pid2nuccore_ids.items()]
    _nuccore_summary = []
    nuccore_summary,failed = edl.esummary(db='nuccore',
                                    ids=nuccore_ids,
                                    result_func=lambda x: Entrez.read(io.StringIO(x)))
    if failed:
        _nuccore_summary,failed = edl.esummary(db='nuccore',
                                ids=failed,
                                result_func=lambda x: Entrez.read(io.StringIO(x)),
                                batch_size=1)
    nid2length = {}
    for result in nuccore_summary+_nuccore_summary:
        nid2length[result['AccessionVersion']] = result['Length'].real
    nid2length = {k:v for k,v in nid2length.items() if k in pid2nuccore_ids.values()}
    
    tqdm.write('%s need to manual adjust....' % len(set(pid2nuccore_ids.values()).difference(set(nid2length))))
    near_end_protein = []
    for idx,row in full_df.iterrows():
        end = row['end']
        start = row['start']
        nuccore_ID = pid2nuccore_ids[idx]
        length = nid2length.get(nuccore_ID,'')
        if length:
            if abs(end - length) <= 100 or start <=100:
                tqdm.write('near end :' + idx)
                near_end_protein.append(idx)
        else:
            return_v = input(f'nuccore ID is {nuccore_ID}; does not get any info, maybe wgs, please check it manually. protein id is {idx}, end at {end}. Y/y for indicated it near the end. ')
            if return_v.lower() == 'y':
                near_end_protein.append(idx)
            else:
                pass
    return near_end_protein
        
def filter_archaea(full_df,remove_nc=False):
    # filter out archaea and NC10
    if remove_nc:
        remained_ids = full_df.index[(full_df.loc[:,'superkingdom']=='Bacteria') & (~full_df.loc[:,'phylum'].str.contains('NC10').fillna(False))]
    else:
        remained_ids = full_df.index[full_df.loc[:,'superkingdom']=='Bacteria']
    return remained_ids


tds = ['nr_retrieve_amoB','nr_retrieve_amoC','with_genome_amoA','nr_retrieve_hao']
for target_dir in tds:
    g = target_dir.split('_')[-1]
    if g in ['amoB','amoC']:
        pro2full_tab = f'{target_dir}/filtered_by_kegg.faa_aln.dir/iqtree.treefile/info_dir/pro2full_info.tab'
        fa = f'{target_dir}/filtered_by_kegg.faa'
    elif g == 'amoA':
        pro2full_tab = f'{target_dir}/used.fasta_aln.dir/iqtree.treefile/info_dir/pro2full_info.tab'
        fa = f'{target_dir}/filtered_by_kegg.faa'
    elif g == 'hao':
        pro2full_tab = f'{target_dir}/filtered_by_kegg.faa_aln.dir/iqtree.no_trim.treefile/info_dir/pro2full_info.tab'
        fa = f'{target_dir}/filtered_by_kegg.faa'
        
    full_df = pd.read_csv(pro2full_tab,sep='\t',index_col=0)
    records = list(SeqIO.parse(fa,format='fasta'))
    near_end_protein = filter_by_relative_pos(full_df)

    remained_B_ids = filter_archaea(full_df,remove_nc=False)
    remained_B_remove_nc_ids = filter_archaea(full_df,remove_nc=True)
    num_ori = len(records)

    remained_records = [_ for _ in records if _.id in remained_B_ids]
    remained_records = [_ for _ in remained_records if _.id not in near_end_protein]
    with open(f'{target_dir}/with_genome_Bacteria_intact.faa','w') as f1:
        SeqIO.write(remained_records,f1,format='fasta-2line')
    print('remained %s fa' % len(remained_records))
    print('original %s fa' % num_ori)

    remained_records = [_ for _ in records if _.id in remained_B_remove_nc_ids]
    remained_records = [_ for _ in remained_records if _.id not in near_end_protein]
    with open(f'{target_dir}/with_genome_Bacteria_drop_NC10_intact.faa','w') as f1:
        SeqIO.write(remained_records,f1,format='fasta-2line')
    print('remained %s fa' % len(remained_records))
    print('original %s fa' % num_ori)

    cmd = f'python3 ~/script/evolution_relative/global_search/build_tree_exe.py {target_dir}/with_genome_Bacteria_intact.faa'
    check_call(cmd,shell=1)

    cmd = f'python3 ~/script/evolution_relative/global_search/build_tree_exe.py {target_dir}/with_genome_Bacteria_drop_NC10_intact.faa'
    check_call(cmd,shell=1)
    ####### after iqtree
    cmd = f'python3 ~/script/evolution_relative/global_search/build_tree_exe.py {target_dir}/with_genome_Bacteria_intact.faa .iqtree.treefile'
    check_call(cmd,shell=1)

    cmd = f'python3 ~/script/evolution_relative/global_search/build_tree_exe.py {target_dir}/with_genome_Bacteria_drop_NC10_intact.faa .iqtree.treefile'
    check_call(cmd,shell=1)
    cmd = f"cp -r {dirname(pro2full_tab)} {target_dir}/with_genome_Bacteria_intact.faa_aln.dir/iqtree.treefile" 
    check_call(cmd,shell=1)
    cmd = f"cp -r {dirname(pro2full_tab)} {target_dir}/with_genome_Bacteria_drop_NC10_intact.faa_aln.dir/iqtree.treefile" 
    check_call(cmd,shell=1)
    cmd = f'python3 ~/script/evolution_relative/global_search/reannotate_tree.py {target_dir}/with_genome_Bacteria_intact.faa_aln.dir/iqtree.treefile {target_dir}/with_genome_Bacteria_drop_NC10_intact.faa_aln.dir/iqtree.treefile'
    check_call(cmd,shell=1)
    
if __name__ == "__main__":
    import sys
    params = sys.argv[1:]
    
    