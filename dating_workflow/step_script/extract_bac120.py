from glob import glob
from subprocess import check_call, check_output
import os
from tqdm import tqdm
from os.path import *
from collections import defaultdict
from Bio import SeqIO
import multiprocessing as mp
from dating_workflow.step_script import _parse_blastp, _parse_hmmscan, _get_tophit
import click


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
        r = list(tqdm(tp.imap(run, params), total=len(params)))


def parse_annotation(odir, top_hit=False, evalue=1e-50):
    # for cdd
    _cdd_match_ids = pfam_ids
    genome2annotate = defaultdict(lambda: defaultdict(list))

    # cdd annotations
    tqdm.write('start to read/parse output files (cdd and tigrfam)')
    cdd_anno_files = glob(join(odir, 'PFAM', '*.out'))
    # tigrfam annotations
    tigrfam_anno_files = glob(join(odir, 'TIGRFAM', '*.out'))
    for ofile in tqdm(tigrfam_anno_files + cdd_anno_files):
        gname = basename(ofile).replace('.out', '')
        locus_dict = _parse_hmmscan(ofile=ofile,
                                    top_hit=top_hit,
                                    filter_evalue=evalue)
        genome2annotate[gname].update(locus_dict)
    return genome2annotate


# extract protein
def write_cog(outdir, genome2cdd, raw_proteins, genome_ids=[], get_type='prot'):
    genome2seq = {}
    if not genome_ids:
        genome_ids = list(genome2cdd)
    gene_ids = set([_ for vl in genome2cdd.values() for _ in vl])
    pdir = dirname(expanduser(raw_proteins))
    if get_type == 'nuc':
        suffix = 'ffn'
    elif get_type == 'prot':
        suffix = 'faa'
    else:
        raise Exception
    if not exists(outdir):
        os.makedirs(outdir)
    tqdm.write('get sequence file')
    for genome_name in tqdm(genome_ids):
        g_dict = genome2cdd[genome_name]
        gfile = f'{pdir}/{genome_name}.faa'
        new_pdir = abspath(dirname(dirname(realpath(gfile))))
        gfile = f"{new_pdir}/tmp/{genome_name}/{genome_name}.{suffix}"

        if exists(gfile):
            _cache = {record.id: record
                      for record in SeqIO.parse(gfile, format='fasta')}
            seq_set = {k: [_cache[_]
                           for _ in v
                           if _ in _cache]
                       for k, v in g_dict.items()}
            genome2seq[genome_name] = seq_set

    # concat/output proteins
    tqdm.write('write out')
    for each_gene in tqdm(gene_ids):
        gene_records = []
        for gname, seq_dict in genome2seq.items():
            get_records = seq_dict.get(each_gene, [])
            for record in get_records:
                record.name = gname
            gene_records += get_records
        unique_cdd_records = []
        [unique_cdd_records.append(record)
         for record in gene_records
         if record.id not in [_.id
                              for _ in unique_cdd_records]]

        with open(join(outdir, f"{each_gene.replace('CDD:','')}.faa"), 'w') as f1:
            SeqIO.write(unique_cdd_records, f1, format='fasta-2line')

# def perform_iqtree(indir):
#     script = expanduser('~/bin/batch_run/batch_mafft.py')
#     run(f"python3 {script} -i {indir} -o {indir}")

#     script = expanduser('~/bin/batch_run/batch_tree.py')
#     run(f"python3 {script} -i {indir} -o {join(dirname(indir),'tree')} -ns newick -use fasttree")


def stats_cog(genome2genes):
    gene_ids = pfam_ids+tigfam_ids

    gene_multi = {g: 0 for g in gene_ids}
    for genome, pdict in genome2genes.items():
        for gene, seqs in pdict.items():
            if len(seqs) >= 2:
                gene_multi[gene] += 1
    gene_Ubiquity = {g: 0 for g in gene_ids}
    for genome, pdict in genome2genes.items():
        for gene, seqs in pdict.items():
            gene_Ubiquity[gene] += 1

    gene2genome_num = {}
    for gene in gene_ids:
        _cache = [k for k, v in genome2genes.items() if v.get(gene, [])]
    # for genome, pdict in genome2genes.items():
        gene2genome_num[gene] = len(_cache)

    return gene_multi, gene_Ubiquity, gene2genome_num


def process_path(path):
    if '~' in path:
        path = expanduser('path')
    if not '/' in path:
        path = './' + path
    path = abspath(path)
    return path


@click.command()
@click.option("-in_p", 'in_proteins',)
@click.option("-in_a", 'in_annotations',)
@click.option("-s", "suffix", default='faa')
@click.option("-o", 'outdir',)
@click.option("-evalue", 'evalue', default=1e-50)
@click.option("-gl", "genome_list", default=None, help="It will read 'selected_genomes.txt', please prepare the file, or indicate the alternative name or path. It could be None. If you provided, you could use it to subset the aln sequences by indicate names.")
def main(in_proteins, suffix, in_annotations, outdir, evalue, genome_list):
    if genome_list is None:
        gids = []
    else:
        gids = open(genome_list).read().split('\n')
        gids = list(set([_ for _ in gids if _]))
    in_proteins = join(in_proteins, '*.'+suffix.strip('.'))
    protein_files = glob(in_proteins)
    gids = []
    if not protein_files:
        exit(f"error input proteins dir {in_proteins}")

    annotate_bac120(protein_files, in_annotations, db_id='tigrfam')
    annotate_bac120(protein_files, in_annotations, db_id='pfam')

    genome2genes = parse_annotation(in_annotations, top_hit=False)
    gene_multi, gene_Ubiquity, gene2genome_num = stats_cog(genome2genes)

    genome2genes = parse_annotation(
        in_annotations, top_hit=True, evalue=evalue)
    write_cog(outdir, genome2genes, in_proteins,
              genome_ids=gids, get_type='prot')


if __name__ == "__main__":
    main()

    # import sys

    # # usage :
    # # extract_cog.py 'raw_genome_proteins/*.faa' ./target_genes ./conserved_protein
    # if len(sys.argv) >= 2:
    #     raw_proteins = process_path(sys.argv[1])
    #     out_cog_dir = process_path(sys.argv[2])
    #     outdir = process_path(sys.argv[3])
    #     protein_files = glob(raw_proteins)
    #     gids = []
    # else:
    #     raw_proteins = expanduser('~/data/nitrification_for/dating_for/raw_genome_proteins/*.faa')
    #     out_cog_dir = expanduser('~/data/nitrification_for/dating_for/bac120_annoate')
    #     outdir = expanduser('~/data/nitrification_for/dating_for/bac120_annoate/seq')
    #     protein_files = glob(raw_proteins)

    # # for tigfam_id in tigfam_ids:
    # annotate_bac120(protein_files, out_cog_dir, db_id='tigrfam')
    # annotate_bac120(protein_files, out_cog_dir, db_id='pfam')

    # genome2genes = parse_annotation(out_cog_dir,top_hit=False)
    # gene_multi,gene_Ubiquity,gene2genome_num = stats_cog(genome2genes)

    # genome2genes = parse_annotation(out_cog_dir,top_hit=True)
    # write_cog(outdir,genome2genes,raw_proteins,genome_ids=gids,get_type='prot')

    # # perform_iqtree(outdir)
