from os.path import *
import os
from subprocess import check_call
from tqdm import tqdm

cdd_odir = join('./genome_protein_files_more','cogRbp_anno')
os.makedirs(cdd_odir,exist_ok=1)
for fa in tqdm(fa_files):
    anno = join(cdd_odir,basename(fa).replace('.faa','.cogrbp'))
    cmd_template = f"blastp -query {fa} -db /home-user/sswang/resource/db/CogRbp/CogRbp -outfmt 6 -num_threads 50 -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 > {anno}" 
    check_call(cmd_template,shell=1)
    