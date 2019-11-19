from glob import glob
from subprocess import check_call
import os
from tqdm import tqdm
from os.path import *
from collections import defaultdict
from Bio import SeqIO
import multiprocessing as mp


cdd_list = "/home-db/pub/protein_db/CDD/cdd.info.dat/cddid_all.tbl"
__file__ = '/home-user/thliao/script/evolution_relative/dating_workflow/extrat_cog.py'
cog_list = join(dirname(__file__),'cog.list')
cog_list = set([_ for _ in open(cog_list).read().split('\n') if _])

cdd_num = []
for row in open(cdd_list,'r'):
    if row.split('\t')[1] in cog_list:
        cdd_num.append("CDD:%s" % row.split('\t')[0])
# ABOVE is the default setting for luolab server.



def run(cmd):
    check_call(cmd,
               shell=True)

def annotate_cog(raw_protein,out_cog_dir):
    params = []
    for f in tqdm(glob(raw_protein)):
        gname = f.split('/')[-1].replace('.faa', '')
        ofile = f'{out_cog_dir}/{gname}.out'
        cmd = f"/home-user/software/blast/latest/bin/rpsblast -query {f} -db /home-db/pub/protein_db/CDD/Cog -max_target_seqs 1 -num_threads 30 -outfmt 6 -evalue 1e-5 -out {ofile}"
        if not os.path.exists(ofile):
            params.append(cmd)
            # check_call(cmd, shell=1)
    with mp.Pool(processes=5) as tp:
        list(tqdm(tp.imap(run,params),total=len(params)))


def extra_cog(cog_out_dir):
    # cog out dir
    genome2cdd = defaultdict(lambda:defaultdict(list))
    for f in tqdm(glob(join(cog_out_dir,'*.out'))):
        genome_name = f.split('/')[-1].replace('.out','')
        for row in open(f,'r'):
            locus = row.split('\t')[0]
            if row.split('\t')[1] in cdd_num:
                genome2cdd[genome_name][row.split('\t')[1]].append(locus)
    return genome2cdd
# tmp (prokka)
# rrna doesn't have protein sequences, pass it
# tmp_list = [expanduser('~/data/nitrification_for/tmp'),
#             expanduser('~/data/nitrification_for/cyano_basal/tmp'),
#             expanduser('~/data/nitrification_for/backbone_others/tmp')]

# for prokka_dir in tmp_list:
#     for tsv in tqdm(glob(join(prokka_dir,'*/*.tsv'))):
#         genome_name = tsv.split('/')[-1].replace('.tsv','')
#         for row in open(tsv):
#             if '16S ribosomal RNA' == row.split('\t')[-1].strip('\n'):
                
#                 genome2cdd[genome_name]['16S'].append(row.split('\t')[0])
#             elif '23S ribosomal RNA' in row.split('\t')[-1].strip('\n'):
#                 genome2cdd[genome_name]['23S'].append(row.split('\t')[0])

# extract protein
def write_cog(outdir,genome2cdd,raw_proteins):
    pdir = dirname(expanduser(raw_proteins))
    for genome_name,pset in tqdm(genome2cdd.items()):
        pfiles = glob(f'{pdir}/{genome_name}.faa')
        if pfiles:
            pfile = pfiles[0]
            _cache = {record.id:record for record in SeqIO.parse(pfile,format='fasta')}
            pset = {k:[_cache[_] for _ in v if _ in _cache] for k,v in pset.items()}
            genome2cdd[genome_name] = pset
        else:
            continue
    # concat/output proteins
    os.makedirs(outdir,exist_ok=1)
    for each_cdd in tqdm(cdd_num):
        cdd_records = []
        for gname, p_d in genome2cdd.items():
            get_records = p_d.get(each_cdd,[])
            if len(get_records) >=1:
                record = get_records[0]
                record.id = gname
                cdd_records.append(record)
        with open(join(outdir,f"{each_cdd.replace('CDD:','')}.faa"),'w') as f1:
            SeqIO.write(cdd_records,f1,format='fasta-2line')
    
              

    from collections import Counter
    for g,v in genome2cdd.items():
        counter_v = Counter(v)
        if [(k,v) for k,v in counter_v.items() if v >=2]:
            print([(k,v) for k,v in counter_v.items() if v >=2])
            
if __name__ == "__main__":
    import sys
    # usage :
    # extract_cog.py 'raw_genome_proteins/*.faa' ./target_genes ./conserved_protein
    if len(sys.argv) >= 2:
        raw_proteins = sys.argv[1]
        out_cog_dir = sys.argv[2]
        outdir = sys.argv[3]
    else:
        raw_proteins = expanduser('~/data/nitrification_for/dating_for/raw_genome_proteins/*.faa')
        out_cog_dir = expanduser('~/data/nitrification_for/dating_for/target_genes')
        outdir = expanduser('~/data/nitrification_for/dating_for/conserved_protein')
    annotate_cog(raw_proteins, out_cog_dir)
    genome2cdd = extra_cog(out_cog_dir)
    write_cog(outdir, genome2cdd,raw_proteins)