import pandas as pd
from os.path import join,exists,basename,dirname,isdir,expanduser
from glob import glob
import os
from subprocess import check_call
from tqdm import tqdm

base_tab = expanduser('~/.cache/ncbi-genome-download/genbank_bacteria_assembly_summary.txt')
indir = './genbank'
odir = './genome_protein_files'
all_g_ids = set([basename(_) for _ in glob(join(indir,'bacteria','*'))])

metadatas = open(base_tab).read().split('\n')
rows = [_ for _ in metadatas if _.split('\t')[0] in all_g_ids]
with open('./metadata.csv','w') as f1:
    f1.write(metadatas[1].strip('# ') + '\n') 
    f1.write('\n'.join(rows))
tmp_dir = './tmp'
prokka_p = "/usr/local/bin/prokka"
if not exists(odir):
    os.makedirs(odir,exist_ok=True)
if not exists(tmp_dir):
    os.makedirs(tmp_dir,exist_ok=True)
    
def run_cmd(cmd):
    check_call(cmd,shell=True,stderr=open('/dev/null','w'),stdout =open('/dev/null','w'))
def get_faa_from_prokka_r(infile,odir,sample_name):
    locustag = sample_name.split('_')[-1].split('.')[0]
    if exists(f'{odir}/{sample_name}/{sample_name}.faa'):
        return f'{odir}/{sample_name}/{sample_name}.faa'
    run_cmd(f'{prokka_p} {infile} --outdir {odir}/{sample_name} --force --prefix {sample_name} --locustag {locustag} --cpus 0 ')
    return f'{odir}/{sample_name}/{sample_name}.faa'


for p_dir in tqdm([_ for _ in glob(join(indir,'**','GCA*'),recursive=True) if isdir(_)]):
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