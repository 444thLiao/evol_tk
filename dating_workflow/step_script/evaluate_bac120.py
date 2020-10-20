

# evaulate faa first

from glob import glob
import os
from os.path import *
from Bio import SeqIO


def check_r(r):
    
    names = r.description.split(' ')
    if len(names)>=2:
        if names[0].split('_')[0] == names[1].split('_')[0] and names[1].split('_')[0]:
            r.description = names[0] + ' ' + ' '.join(names[2:])
            return r
    else:
        return None



faa_files = glob('./modified_data/genome_protein_files/*.faa')
for faa in tqdm(faa_files):
    records = list(SeqIO.parse(faa,format='fasta'))
    new_records = []
    contains_wrong = False
    for r in records:
        if check_r(r) is not None:
            contains_wrong = True
        new_records.append(r)


faa_files = glob("./modified_data/genome_protein_files/*.faa")
for faa in faa_files:
    if getsize(faa)==0:
        os.system(f'rm {faa}')



for faa in tqdm(glob('./modified_data/prokka_o/*/*.faa')):
    target_dir =  "./modified_data/genome_protein_files/"
    if exists(faa.replace('.faa','.gbk')) and not exists(join(target_dir,basename(faa))) and getsize(faa)!=0:
        os.system(f"ln -sf `realpath {faa}`  {target_dir}")
    #records = list(SeqIO.parse(faa,format='fasta'))
rm_ids = []
for faa in glob('./modified_data/prokka_o/*/*.faa'):
    target_dir =  "./modified_data/direct_protein_files/"
    if not exists(faa.replace('.faa','.gbk')):
        rm_ids.append(dirname(faa))
    
    with open(f'./modified_data/direct_protein_files/{basename(faa)}','w') as f1:
        SeqIO.write(records,f1,format='fasta-2line')


# format it
from glob import glob
from os.path import *
import os
indir = './genome_protein_files'
for faa in glob(join(indir,'*.faa')):
    if islink(faa):
        os.system(f'rm {faa}; cp `realpath {faa}` {indir}/')

for aid in fna_ids:
    f ="/mnt/maple/thliao/data/NCBI/modified_data/genome_protein_files"
    if not exists(f'./genome_protein_files/{aid}.faa'):
        os.system(f'cp {f}/{aid}.faa {indir}/')