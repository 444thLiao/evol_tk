
from api_tools.tk import add_gene,os,SeqIO


for g in comp2name:
    gbk = f"/mnt/ivy/thliao/project/coral_ruegeria/data_processing/pub_dataset/prokka_o/{g}/{g}.gbk"
    print(gbk)
    allcontig = add_gene(gbk)
    os.system(f"mv {gbk} {gbk.replace('.gbk','.old.gbk')}")
    with open(gbk,'w') as f1:
        SeqIO.write(allcontig,f1,'genbank')    
        

from glob import glob
from tqdm import tqdm

gid2fna = {_.split('/')[-3]:_ 
           for _ in glob('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/*/09_prokka/*.fna')
           if '_' not in _.split('/')[-1] }
from Bio.SeqFeature import SeqFeature, FeatureLocation
for g in tqdm(gid2fna):
    gbk = f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/{g}/09_prokka/{g}.gbk"
    allcontig = add_gene(gbk)
    os.system(f"mv {gbk} {gbk.replace('.gbk','.wrong.gbk')}")
    with open(gbk,'w') as f1:
        SeqIO.write(allcontig,f1,'genbank')
        