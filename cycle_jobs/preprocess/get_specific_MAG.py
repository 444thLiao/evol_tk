import pandas as pd
from os.path import join,dirname,exists,basename
from glob import glob
from subprocess import check_call
import os
from tqdm import tqdm
target_file = '/mnt/home-backup/thliao/metagenomes/update_0928_nitrification/confirmed_locus2info.tsv'
genome_dir = '/mnt/home-backup/thliao/metagenomes/concat_all/prokka_o/'
#odir = '/home-user/thliao/project/nitrogen_cycle/nitrification/reference_genomes/genome_protein_files'
#odir = '/home-user/thliao/project/nitrogen_cycle/nitrification/reference_genomes/genome_protein_files_more'

odir = '/home-user/thliao/project/nitrogen_cycle/nitrification/reference_genomes/test/fna'
if not exists(odir):
    os.makedirs(odir,)

name2ko = {'amoA': 'K10944',
           'amoB': 'K10945',
           'amoC': 'K10946',
           'hao': 'K10535',
           'nxrA': 'K00370',
           'nxrB': 'K00371'}
info_df = pd.read_csv(target_file,sep='\t')
sub_df = info_df.loc[info_df.loc[:,'Gene name(N metabolism)'].isin(name2ko),:]
sample_names = sub_df.loc[:,'sample name']
for sname in tqdm(sample_names):
    p_dir = join(genome_dir,sname)
    faa = [_ for _ in glob(join(p_dir,'*.fna'))][0]
    if not exists(join(odir,basename(faa))):
        check_call(f"ln -s {faa} {odir}",shell=True)
        

