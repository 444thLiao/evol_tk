from subprocess import check_call
from glob import glob 
from os.path import join,dirname,exists,basename
import os
from tqdm import tqdm
from collections import defaultdict


indir = './genome_protein_files'
odir = './annotated'

kofam_scan = '/home-user/thliao/software/kofamscan/exec_annotation'

def run(cmd):
    check_call(cmd,shell=True)

for infile in tqdm(glob(join(indir,'*.faa'))):
    sname = basename(infile).replace('.faa','')
    ofile = join(odir,sname+'.out')
    if not exists(odir):
        os.makedirs(odir,exist_ok=True)
    if not exists(ofile):
        run(f"{kofam_scan} -o {ofile} --cpu 64 -f mapper-one-line --no-report-unannotated {infile}")
