from glob import glob
from subprocess import check_call, check_output
import os
from tqdm import tqdm
from os.path import *
from collections import defaultdict
from Bio import SeqIO
import multiprocessing as mp

pfam_db = '/home-user/thliao/data/protein_db/bac120/Pfam.v32.sub6.hmm'
tigfam_db = '/home-user/thliao/data/protein_db/bac120/TIGRFAMv14_sub114.hmm'

__file__ = '/home-user/thliao/script/evolution_relative/dating_workflow/step_script/extrat_bac120.py'
bac120_list = join(dirname(__file__), 'bac120.tsv')
id_list = [row.split('\t')[0] for row in open(bac120_list) if row]
id_list = id_list[1:]
pfam_ids = [_ for _ in id_list if _.startswith('PF0')]
tigfam_ids = [_ for _ in id_list if _.startswith('TIGR')]


# ABOVE is the default setting for luolab server.


def run(cmd):
    check_call(cmd,
               shell=True,
               stdout=open('/dev/null', 'w'),
               )


def annotate_bac120(protein_files, odir, db_id='pfam'):
    params = []
    if not exists(odir):
        os.makedirs(odir)
    for pfile in protein_files:
        gname = basename(pfile).replace('.faa', '')
        if db_id == 'pfam':
            ofile = f'{odir}/PFAM/{gname}.out'
            cmd = f"/usr/local/bin/hmmscan --tblout {ofile} --acc --noali --notextw --cpu 40 {pfam_db} {pfile}"
        elif db_id == 'tigrfam':
            ofile = f'{odir}/TIGRFAM/{gname}.out'
            cmd = f"/usr/local/bin/hmmscan --tblout {ofile} --acc --noali --notextw --cpu 40 {tigfam_db} {pfile}"
        # else:
        #     ofile = f'{odir}/{db_id}/{gname}.out'
        #     assert exists(f"{tigfam_db_dir}/{db_id}.HMM")
        #     cmd = f"/usr/local/bin/hmmscan --tblout {ofile} --acc --noali --notextw --cpu 40 {tigfam_db_dir}/{db_id}.HMM {pfile}"
        if not exists(ofile):
            if not exists(dirname(ofile)):
                os.makedirs(dirname(ofile))
            params.append(cmd)
            # check_call(cmd, shell=1)
    # print(params)
    with mp.Pool(processes=5) as tp:
        r = list(tqdm(tp.imap(run, params),total=len(params)))


def parse_grep_result(ofile,filter_evalue,get_top=False):
    # for grep against hmm scan
    gid2locus = defaultdict(list)

    for row in open(ofile, 'r'):
        if row.startswith('#'):
            continue
        r = row.split(' ')
        r = [_ for _ in r if _]

        gene_id = r[1]
        locus_tag = r[2]
        evalue = float(r[4])
        if filter_evalue and evalue <= filter_evalue:
            if not get_top:
                gid2locus[gene_id].append(locus_tag)
            else:
                gid2locus[gene_id].append((locus_tag, evalue))
        else:
            gid2locus[gene_id].append((locus_tag, evalue))
    # if get_top:
    #     gid2locus = {gene_id:[list(sorted(v_list,
    #                                       key=lambda x:x[1]))[0][0]]
    #                  for gene_id,v_list in gid2locus.items()}
    return gid2locus


def extra_genes(odir, evalue=1e-10, mode='multiple'):
    genome2gene_id = defaultdict(dict)
    odir_pfam = join(odir, 'PFAM')
    odir_tigrfam = join(odir, 'TIGRFAM')
    outfiles = glob(join(odir_pfam, '*.out')) + glob(join(odir_tigrfam, '*.out'))
    tqdm.write('start to read/parse output files')
    for outfile in tqdm(outfiles):
        gname = basename(outfile).replace('.out', '')
        if mode == 'multiple':
            gid2locus = parse_grep_result(outfile, filter_evalue=evalue)
        elif mode == 'top':
            gid2locus = parse_grep_result(outfile, filter_evalue=None,get_top=True)
        # else:
        #     return
        genome2gene_id[gname].update(gid2locus)
    if mode == 'top':
        for gname in genome2gene_id.keys():
            gid2locus = genome2gene_id[gname]
            gid2locus = {gene_id: [list(sorted(v_list,
                                               key=lambda x: x[1]))[0][0]]
                         for gene_id, v_list in gid2locus.items()}
            genome2gene_id[gname] = gid2locus
    return genome2gene_id


# extract protein
def write_genes_multiple(outdir, genome2gene_id, protein_files):
    pdir = dirname(expanduser(protein_files[0]))
    genome2gene_seq = {}
    for genome_name, pdict in tqdm(genome2gene_id.items()):
        pfiles = glob(f'{pdir}/{genome_name}.faa')
        if pfiles:
            pfile = pfiles[0]
            _cache = {record.id: record
                      for record in SeqIO.parse(pfile, format='fasta')}
            pset = {k: [_cache[_]
                        for _ in v
                        if _ in _cache]
                    for k, v in pdict.items()}
            genome2gene_seq[genome_name] = pset
        else:
            continue
    all_ids = set([_ for v in genome2gene_seq.values() for _ in v.keys()])
    # concat/output proteins
    if not exists(outdir):
        os.makedirs(outdir)
    for gene_id in tqdm(all_ids):
        gene_records = []
        for gname, p_d in genome2gene_seq.items():
            get_records = p_d.get(gene_id, [])
            if len(get_records) >= 1:
                # record = get_records
                for record in get_records:
                    record.name = gname
                gene_records += get_records
        with open(join(outdir, f"{gene_id}.faa"), 'w') as f1:
            SeqIO.write(gene_records, f1, format='fasta-2line')
    return genome2gene_seq

def perform_iqtree(indir):
    script = expanduser('~/bin/batch_run/batch_mafft.py')
    run(f"python3 {script} -i {indir} -o {indir}")

    script = expanduser('~/bin/batch_run/batch_tree.py')
    run(f"python3 {script} -i {indir} -o {join(dirname(indir),'tree')} -ns newick -use fasttree")


def stats_cog(genome2genes):
    gene_multi = defaultdict(int)
    for genome, pdict in genome2genes.items():
        for gene, seqs in pdict.items():
            if len(seqs) >= 2:
                gene_multi[gene] += 1
    gene_Ubiquity = defaultdict(int)
    for genome, pdict in genome2genes.items():
        for gene, seqs in pdict.items():
            gene_Ubiquity[gene] += 1
    return gene_multi,gene_Ubiquity

def process_path(path):
    if '~' in path:
        path = expanduser('path')
    if not '/' in path:
        path = './' + path
    path = abspath(path)
    return path

if __name__ == "__main__":
    import sys

    # usage :
    # extract_cog.py 'raw_genome_proteins/*.faa' ./target_genes ./conserved_protein
    if len(sys.argv) >= 2:
        raw_proteins = process_path(sys.argv[1])
        out_cog_dir = process_path(sys.argv[2])
        outdir = process_path(sys.argv[3])
        protein_files = glob(raw_proteins)
    else:
        raw_proteins = expanduser('~/data/nitrification_for/dating_for/raw_genome_proteins/*.faa')
        out_cog_dir = expanduser('~/data/nitrification_for/dating_for/bac120_annoate')
        outdir = expanduser('~/data/nitrification_for/dating_for/bac120_annoate/seq')

        protein_files = glob(raw_proteins)
    # for tigfam_id in tigfam_ids:
    annotate_bac120(protein_files, out_cog_dir, db_id='tigrfam')
    annotate_bac120(protein_files, out_cog_dir, db_id='pfam')

    genome2genes = extra_genes(out_cog_dir,mode='multiple')
    gene_multi,gene_Ubiquity = stats_cog(genome2genes)

    genome2genes = extra_genes(out_cog_dir, mode='top')
    genome2gene_seq = write_genes_multiple(outdir, genome2genes,protein_files)

    perform_iqtree(outdir)
