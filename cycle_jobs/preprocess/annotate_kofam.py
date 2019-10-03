# Only annoate reference genome for get standard genes

from subprocess import check_call
from glob import glob 
from os.path import join,dirname,exists,basename
import os
from tqdm import tqdm
from collections import defaultdict
import pandas as pd


genome_info = './genome_info_full.xlsx'
g_df = pd.read_excel(genome_info, index_col=0)
g_df = g_df.loc[g_df.loc[:,'used']!='no',:]
all_g_ids = list(g_df.index)

indir = './genome_protein_files'
odir = './annotated'

kofam_scan = '/home-user/thliao/software/kofamscan/exec_annotation'

def run(cmd):
    check_call(cmd,shell=True)

for infile in tqdm(glob(join(indir,'*.faa'))):
    sname = basename(infile).replace('.faa','')
    ofile = join(odir,sname+'.out')
    if sname in all_g_ids:
        if not exists(odir):
            os.makedirs(odir,exist_ok=True)
        if not exists(ofile):
            run(f"{kofam_scan} -o {ofile} --cpu 64 -f mapper-one-line --no-report-unannotated {infile}")
