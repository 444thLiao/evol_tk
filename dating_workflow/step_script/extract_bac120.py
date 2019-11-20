from glob import glob
from subprocess import check_call
import os
from tqdm import tqdm
from os.path import *
from collections import defaultdict
from Bio import SeqIO
import multiprocessing as mp

pfam_db = '/home-db/pub/protein_db/Pfam/Pfam.v32.0/Pfam-A.hmm'
tigfam_db_dir = '/home-user/thliao/data/protein_db/TIGRFAM/'


__file__ = '/home-user/thliao/script/evolution_relative/dating_workflow/step_script/extrat_bac120.py'
bac120_list = join(dirname(__file__), 'cog.list')
id_list = [row.split('\t')[0] for row in open(bac120_list) if _]
id_list = id_list[1:]
pfam_ids = [_ for _ in id_list if _.startswith('PF0')]
tigfam_ids = [_ for _ in id_list if _.startswith('TIGR')]
# ABOVE is the default setting for luolab server.


def run(cmd):
    check_call(cmd,
               shell=True)


def annotate_gene(raw_protein, out_cog_dir,type='pfam'):
    params = []
    for f in tqdm(glob(raw_protein)):
        gname = f.split('/')[-1].replace('.faa', '')
        if type == 'pfam':
            ofile = f'{out_cog_dir}/{gname}.out'
            cmd = f"hmmscan --tblout {ofile} --acc --noali --notextw --cpu 40 {pfam_db}  {f}"
        else:
            ofile = f'{out_cog_dir}/{gname}.{type}.out'
            assert exists(f"{tigfam_db_dir}/{type}.HMM")
            cmd = f"hmmscan --tblout {ofile} --acc --noali --notextw --cpu 40 {tigfam_db_dir}/{type}.HMM  {f}"
        if not os.path.exists(ofile):
            params.append(cmd)
            # check_call(cmd, shell=1)
    with mp.Pool(processes=5) as tp:
        list(tqdm(tp.imap(run, params), total=len(params)))


def extra_gene(cog_out_dir):
    # cog out dir
    genome2cdd = defaultdict(lambda: defaultdict(list))
    for f in tqdm(glob(join(cog_out_dir, '*.out'))):
        genome_name = f.split('/')[-1].replace('.out', '')
        for row in open(f, 'r'):
            locus = row.split('\t')[0]
            if row.split('\t')[1] in cdd_num:
                genome2cdd[genome_name][row.split('\t')[1]].append(locus)
    return genome2cdd


# extract protein
def write_cog_multiple(outdir, genome2cdd, raw_proteins):
    pdir = dirname(expanduser(raw_proteins))
    for genome_name, pdict in tqdm(genome2cdd.items()):
        pfiles = glob(f'{pdir}/{genome_name}.faa')
        if pfiles:
            pfile = pfiles[0]
            _cache = {record.id: record
                      for record in SeqIO.parse(pfile, format='fasta')}
            pset = {k: [_cache[_]
                        for _ in v
                        if _ in _cache]
                    for k, v in pdict.items()}
            genome2cdd[genome_name] = pset
        else:
            continue
    # concat/output proteins
    if not exists(outdir):
        os.makedirs(outdir)
    for each_cdd in tqdm(cdd_num):
        cdd_records = []
        for gname, p_d in genome2cdd.items():
            get_records = p_d.get(each_cdd, [])
            if len(get_records) >= 1:
                # record = get_records
                for record in get_records:
                    record.name = gname
                cdd_records += get_records
        with open(join(outdir, f"{each_cdd.replace('CDD:', '')}.faa"), 'w') as f1:
            SeqIO.write(cdd_records, f1, format='fasta-2line')


def perform_iqtree(outdir):
    script = expanduser('~/bin/batch_run/batch_mafft.py')
    run(f"python3 {script} -i {outdir} -o {outdir}")

    script = expanduser('~/bin/batch_run/batch_tree.py')
    run(f"python3 {script} -i {outdir} -o {outdir} -ns newick -use fasttree")


def stats_cog(genome2cdd, outdir):
    cog_multi = defaultdict(list)
    for g, _d in genome2cdd.items():
        for cog, v in _d.items():
            if len(v) >= 2:
                cog_multi[cog].append(g)


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
    write_cog_multiple(outdir, genome2cdd)
    perform_iqtree(outdir)