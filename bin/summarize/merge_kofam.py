"""
This script is mainly for merge output by all versus each KO
It is suitble for hmm output which query is KO fam, and the subject is your sequence.

"""

from tqdm import tqdm
from glob import glob
from collections import defaultdict
import pandas as pd
import click
from os.path import *
import os
def advanced_dict(in_list):
    used_locus = {}
    _cache = {}
    for key,v,evalue in in_list:
        if key in _cache:
            if v not in used_locus:
                _cache[key] += f',{v}'
            else:
                if evalue <= used_locus[v]:
                    used_locus[v] = evalue
                    _cache[key] += f',{v}'
        _cache[key] = v
    return _cache

@click.command()
@click.option("-i",'tab_dir')
@click.option("-o",'output_dir')
@click.option("-e",'evalue',default=1e-20,help="filter out evalue below it [default:1e-20]")
@click.option("-no_c","not_convert")
def main(tab_dir,evalue,output_dir,not_convert):
    gid2locus2ko = defaultdict(list)
    tab_files = glob(join(tab_dir,'*.tab'))
    if not tab_files:
        exit(f"wrong -i, no {join(tab_dir,'*.tab')} exits")
    for hf in tqdm(tab_files):
        for row in open(hf):
            if row.startswith('#'):
                continue
            r = row.split(' ')
            r = [_ for _ in r if _]
            gene_id = r[0]
            ko = r[2]
            evalue = float(r[4])
            gid2locus2ko[convert_genome_ID_rev(gene_id)].append((gene_id,ko,evalue))
    
    gid2locus2ko = {k:advanced_dict([(_[1],_[0],_[2]) 
                            for _ in v if _[2]<=evalue]) 
                    for k,v in tqdm(gid2locus2ko.items()) }
    tmp_df = pd.DataFrame.from_dict(gid2locus2ko,orient='index')
    if not exists(output_dir):
        os.makedirs(output_dir)
    tmp_df.to_csv(join(output_dir,'merged_hmm_info.tab'),sep='\t',index=1)
    
    
if __name__ == "__main__":
    main()