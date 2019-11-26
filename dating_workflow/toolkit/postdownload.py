"""
It mainly for running prokka after download genomes.
It includes
1. check the missing ID
2. perform prokka and rename the ID
3. generate metadata.csv
"""
import pandas as pd
from os.path import join,exists,basename,dirname,isdir,expanduser
from glob import glob
import os
from subprocess import check_call
from tqdm import tqdm
import sys

def convert_genome_ID(genome_ID):
    # for GCA_900078535.2
    # it will return 
    return genome_ID.split('_')[-1].replace('.','v')


def run_cmd(cmd):
    check_call(cmd,shell=True,
               stderr=open('/dev/null','w'),
               stdout =open('/dev/null','w'))
    
def get_faa_from_prokka_r(infile,odir,sample_name,return_cmd=False):
    locustag = convert_genome_ID(sample_name)
    oprefix = f"{odir}/{sample_name}/{sample_name}"
    if exists(f'{oprefix}.faa'):
        return f'{oprefix}.faa'
    cmd = f'{prokka_p} {infile} --outdir {odir}/{sample_name} --force --prefix {sample_name} --locustag {locustag} --cpus 0 '
    if return_cmd:
        return cmd
    run_cmd(cmd)
    return f'{oprefix}.faa'

def cli(indir,odir=None):
    if odir is None:
        odir = './genome_protein_files'
    all_dir = [_ for _ in glob(join(indir,'**','GCA*'),recursive=True) 
               if isdir(_)]
    
    for p_dir in tqdm(all_dir):
        p_files = glob(join(p_dir,'*.faa.gz'))
        if not p_files:
            fna_file = glob(join(p_dir,'*.fna.gz'))[0]
            new_fna = fna_file.replace('.gz','')
            if not exists(new_fna):
                run_cmd(f'gunzip -d -c {fna_file} > {new_fna}')
            p_files = [get_faa_from_prokka_r(new_fna,tmp_dir,basename(dirname(fna_file)))]
        p_file = p_files[0]
        ofile = join(odir,basename(p_dir)) +'.faa'
        if (not exists(ofile)) and p_file.endswith('.gz') and exists(p_file):
            run_cmd(f'gunzip -d -c {p_file} >{ofile}')
        elif (not exists(ofile)) and exists(p_file):
            run_cmd(f'cat {p_file} >{ofile}')
        else:
            #print(p_file,ofile)
            pass
        
if __name__ == "__main__":
    if len(sys.argv) >= 2:
        indir = sys.argv[1]
        odir = sys.argv[2] #'./genome_protein_files'
    else:
        indir = './genbank'
        odir = './genome_protein_files'
    default_name = './assembly_ids.list'
    base_tab = expanduser('~/.cache/ncbi-genome-download/genbank_bacteria_assembly_summary.txt')

    all_g_ids = set([basename(_) 
                     for _ in glob(join(indir,'bacteria','*'))])
    # from downloaded dir
    all_ids = open(default_name).read().split('\n')
    all_ids = [_ for _ in all_ids if _]
    all_ids = set(all_ids)
    # from id list

    metadatas = open(base_tab).read().split('\n')
    rows = [_ 
            for _ in metadatas 
            if _.split('\t')[0] in all_g_ids]
    with open('./metadata.csv','w') as f1:
        f1.write(metadatas[1].strip('# ') + '\n') 
        f1.write('\n'.join(rows))
    if set(all_ids) != set(all_g_ids):
        print('different id, missing ids: ' + '\n'.join(all_ids.difference(all_g_ids)))


    tmp_dir = './tmp'
    prokka_p = "/usr/local/bin/prokka"
    if not exists(odir):
        os.makedirs(odir,exist_ok=True)
    if not exists(tmp_dir):
        os.makedirs(tmp_dir,exist_ok=True)
    cli(indir,odir)