import pandas as pd
from Bio import SeqIO
import copy
from os.path import exists,join
import os
from glob import glob
from bin.multiple_sbatch import sbatch_all
import time
from subprocess import check_call
# ! build phylogeny based on orthofinder results (removing most of the sister groups)

def get_seqs(og_df, odir, og_list,locus2seq):
    if not exists(odir):
        os.system(f"mkdir -p {odir}")
    og2genome2locus = og_df[og_list].to_dict()
    for og, _d in og2genome2locus.items():
        seqs = []
        for genome, locus in _d.items():
            if not pd.isna(locus):
                r = locus2seq[locus]
                r = copy.deepcopy(r)
                r.id = genome
                r.name = r.description = ''
                seqs.append(r)
        with open(join(odir, f"{og}.faa"), 'w') as f1:
            SeqIO.write(seqs, f1, 'fasta-2line')
            
def main(og_path,odir=None):
    og_df = pd.read_csv(og_path,sep='\t',index_col=0,low_memory=False)
    og_df = og_df.T
    bin_df = og_df.applymap(lambda x: 0 if pd.isna(x) else len(x.split(',')))
    single_copy_all_og = og_df.columns[(bin_df == 1).all()]
    #print(len(single_copy_all_og))
    # not_dup_og = og_df.columns[~(bin_df > 1).any()]
    #print(len(not_dup_og))
    # 17381
    # not_dup_bin_df = bin_df.loc[:, not_dup_og]
    # single_copy_100_og = not_dup_bin_df.columns[(
    #     not_dup_bin_df.sum() >= not_dup_bin_df.shape[0]*1)]
    #print(len(single_copy_100_og))
    all_ids = list(og_df.index)
    locus2seq = {}
    for genome in all_ids:
        locus2seq.update({r.id: r
                        for r in SeqIO.parse(f"{og_path.rsplit('/',4)[0]}/infaa/{genome}.faa", 'fasta')})
    if odir is None:
        odir = f"{og_path.rsplit('/',4)[0]}/OGtree/{len(all_ids)}genomes"
    get_seqs(og_df, odir, single_copy_all_og, locus2seq)
    with open(join(odir, 'genome.list'), 'w') as f1:
        f1.write('\n'.join([f'{g}\t{g}' for g in all_ids]))

    cmds = []
    for faa in glob(join(odir, '*.faa')):
        aln = faa.replace('.faa','.aln')
        cmd = f"mafft --anysymbol --auto {faa} > {aln}; trimal -in {aln} -out {aln.replace('.aln','.trimal')} -automated1"
        cmds.append(cmd)
    sbatch_all(cmds, batch_size=10, prefix_name='mafft')

    while 1:
        if len(glob(join(odir, '*.trimal')))==len(cmds):
            cmd = f"python3 /home-user/thliao/script/evol_tk/dating_workflow/bin/concat_aln.py -i {odir} -o  {odir}/concat/concat.trimal -s trimal -gl {join(odir,'genome.list')} -ct partition -simple"
            check_call(cmd, shell=1)
            break
        time.sleep(60)
    
    if not exists(f"{odir}/iqtree/"):
        os.system(f"mkdir -p {odir}/iqtree")
    cmd = f"iqtree -nt 20 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre {odir}/iqtree/OG_concat -s {odir}/concat/concat.trimal"
    sbatch_all([cmd], thread_per_tasks=20,prefix_name='OG_all', fixed_cluster='others')


## removing OG that are identified as pseudogenes.
# g2func_df = {}
# for gid in tqdm(gid2nano_fna):
#     _df = pd.read_excel(f'/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/summarized_tables/{gid}_{g2pop[gid]}_functional_genes.xlsx')
#     g2func_df[gid] = _df.loc[~_df['locus tag'].isna(),:]
# locus2pseudo = set()
# for gid, _df in g2func_df.items():
#     _locus = list(_df.loc[_df['pseudogenized']!='NOT','locus tag'].values)    
#     locus2pseudo = locus2pseudo.union(set(_locus))
   
# target_g = [_ for _ in g2pop if g2pop[_] in ['MC10','MC46']]
# sub_bin_df = bin_df.loc[target_g,:]
# single_copy_all_og = og_df.columns[(sub_bin_df == 1).all()]
# print(len(single_copy_all_og))
# sub_og_df = og_df.loc[target_g,single_copy_all_og]
# newbin_df = sub_og_df.applymap(lambda x: 1 if x in locus2pseudo else 0)
# og_list = newbin_df.columns[newbin_df.sum()==0]
# odir = '/mnt/ivy/thliao/project/coral_ruegeria/data_processing/phylogeny/OG_phylogeny/MC10_46_woPseudo_1746OGs/'

################
