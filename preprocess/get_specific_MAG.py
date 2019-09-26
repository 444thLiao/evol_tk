import pandas as pd
from os.path import join,dirname,exists,basename
from glob import glob
from subprocess import check_call

target_file = './confirmed_locus2info.tsv'
genome_dir = '/mnt/home-backup/thliao/metagenomes/concat_all/prokka_o/'
odir = '/home-user/thliao/project/nitrogen_cycle/nitrification/reference_genomes/genome_protein_files'
query_g = 'mo'  # including all pmo and amo


info_df = pd.read_csv(target_file,sep='\t')
sub_df = info_df.loc[info_df.loc[:,'Gene name(N metabolism)'].str.lower().str.contains('mo'),:]
sample_names = sub_df.loc[:,'sample name']
for sname in sample_names:
    p_dir = join(genome_dir,sname)
    faa = [_ for _ in glob(join(p_dir,'*.faa'))][0]
    if not exists(join(odir,basename(faa))):
        check_call(f"ln -s {faa} {odir}",shell=True)