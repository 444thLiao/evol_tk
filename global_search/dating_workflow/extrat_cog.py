from glob import glob
from subprocess import check_call
import os
from tqdm import tqdm
from os.path import *
from collections import defaultdict
cdd_list = "/home-db/pub/protein_db/CDD/cdd.info.dat/cddid_all.tbl"
cog_list = join(dirname(__file__),'cog.list')
cog_list = set([_ for _ in open(cog_list).read().split('\n') if _])

cdd_num = []
for row in open(cdd_list,'r'):
    if row.split('\t')[1] in cog_list:
        cdd_num.append("CDD:%s" % row.split('\t')[0])

# cog out dir
cog_out_dir = expanduser('~/data/nitrification_for/dating/target_genes')
genome2cdd = defaultdict(lambda:defaultdict(list))
for f in tqdm(glob(join(cog_out_dir,'*.out'))):
    genome_name = f.split('/')[-1].replace('.out','')
    for row in open(f,'r'):
        locus = row.split('\t')[0]
        if row.split('\t')[1] in cdd_num:
            genome2cdd[genome_name][row.split('\t')[1]].append(locus)

# tmp (prokka)
tmp_list = [expanduser('~/data/nitrification_for/tmp'),
            expanduser('~/data/nitrification_for/cyano_basal/tmp'),
            expanduser('~/data/nitrification_for/backbone_others/tmp')]

for prokka_dir in tmp_list:
    for tsv in tqdm(glob(join(prokka_dir,'*/*.tsv'))):
        genome_name = tsv.split('/')[-1].replace('.tsv','')
        for row in open(tsv):
            if '16S ribosomal RNA' == row.split('\t')[-1].strip('\n'):
                
                genome2cdd[genome_name]['16S'].append(row.split('\t')[0])
            elif '23S ribosomal RNA' in row.split('\t')[-1].strip('\n'):
                genome2cdd[genome_name]['23S'].append(row.split('\t')[0])
                
                
from collections import Counter
for g,v in genome2cdd.items():
    counter_v = Counter(v)
    if [(k,v) for k,v in counter_v.items() if v >=2]:
        print([(k,v) for k,v in counter_v.items() if v >=2])
            
if __name__ == "__main__":
    
    for f in tqdm(glob('../raw_genome_proteins/*.faa')):
        gname = f.split('/')[-1].replace('.faa','')
        ofile = f'./{gname}.out'
        cmd = f"/home-user/software/blast/latest/bin/rpsblast -query {f} -db /home-db/pub/protein_db/CDD/Cog -max_target_seqs 1 -num_threads 30 -outfmt 6 -evalue 1e-5 -out ./{gname}.out"
        if not os.path.exists(ofile):
            check_call(cmd,shell=1)