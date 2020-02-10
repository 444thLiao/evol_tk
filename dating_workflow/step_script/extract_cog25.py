"""
extract 25 proteins for dating analysis
"""
import multiprocessing as mp
import os
from collections import defaultdict
from glob import glob
from os.path import *
from subprocess import check_call

import click
from Bio import SeqIO
from tqdm import tqdm

from dating_workflow.step_script import _parse_blastp

resource_dir = "/home-user/thliao/data/protein_db/dating_resource"
cog_db = f"{resource_dir}/cog25_rps/sing"
cdd_tbl = f"{resource_dir}/cog/cddid_all.tbl"
list27_genes = f"{resource_dir}/single.cog.list"
full_text = open(list27_genes).read().split('\n')
cog_list = set([_.split('\t')[0]
                for _ in full_text
                if _])

num_cdd2name = {}
cdd_num = defaultdict(list)
for row in open(cdd_tbl, 'r'):
    rows = row.split('\t')
    if rows[1] in cog_list:
        cdd_num[rows[1]].append("CDD:%s" % rows[0])
        num_cdd2name[rows[0]] = rows[2]
cdd_num.pop('TIGR00487')


# TIGRFAM_db = f"{resource_dir}/TIGRFAM_v14/TIGR00487.HMM"
# ABOVE is the default setting for luolab server.

def run(cmd):
    check_call(cmd,
               shell=True,
               stdout=open('/dev/null', 'w'),
               stderr=open('/dev/null', 'w'))


def annotate_cog(protein_file_list, cog_out_dir):
    params = []
    for f in protein_file_list:
        gname = basename(f).replace('.faa', '')
        # for cdd
        ofile = f'{cog_out_dir}/{gname}.out'
        cmd = f"/home-user/software/blast/latest/bin/rpsblast -query {f} -db {cog_db} -max_target_seqs 1 -num_threads 10 -outfmt 6 -evalue 1e-3  -out {ofile}"
        if not os.path.exists(ofile):
            if not exists(dirname(ofile)):
                os.makedirs(dirname(ofile))
            params.append(cmd)
        # for tigrfam
        # ofile = f'{cog_out_dir}/TIGRFAM/{gname}.out'
        # cmd = f"hmmscan --tblout {ofile} --acc --noali --notextw --cpu 10 {TIGRFAM_db} {f}"
        # if not os.path.exists(ofile):
        #     if not exists(dirname(ofile)):
        #         os.makedirs(dirname(ofile))
        #     params.append(cmd)
    with mp.Pool(processes=5) as tp:
        list(tqdm(tp.imap(run, params), total=len(params)))


def parse_annotation(cog_out_dir, top_hit=False, evalue=1e-3):
    # for cdd
    # _cdd_match_ids = set([_ for vl in cdd_num.values() for _ in vl])
    genome2cdd = defaultdict(lambda: defaultdict(list))

    # cdd annotations
    tqdm.write('start to read/parse output files')
    cdd_anno_files = glob(join(cog_out_dir, '*.out'))
    for ofile in tqdm(cdd_anno_files):
        gname = basename(ofile).replace('.out', '')
        locus_dict = _parse_blastp(ofile=ofile,
                                   match_ids=[],
                                   top_hit=top_hit,
                                   filter_evalue=evalue)
        genome2cdd[gname].update(locus_dict)
    # tigrfam annotations
    # tigrfam_anno_files = glob(join(cog_out_dir,'TIGRFAM','*.out'))
    # for ofile in tqdm(tigrfam_anno_files):
    #     gname = basename(ofile).replace('.out','')
    #     locus_dict = _parse_hmmscan(ofile=ofile,
    #                                top_hit=top_hit)
    #     genome2cdd[gname].update(locus_dict)
    return genome2cdd


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
        if suffix == 'faa':
            # important bugs!!!!! fixed
            new_gfile = gfile
        else:
            new_gfile = f"{new_pdir}/tmp/{genome_name}/{genome_name}.{suffix}"

        if exists(new_gfile):
            _cache = {record.id: record
                      for record in SeqIO.parse(new_gfile, format='fasta')}
            seq_set = {k: [_cache[_]
                           for _ in v
                           if _ in _cache]
                       for k, v in g_dict.items()}
            genome2seq[genome_name] = seq_set
        else:
            # not with prokka annotations
            print('not annotated with prokka')
            if not gfile.endswith(suffix):
                print(f'not {suffix},past it')
                continue
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

        with open(join(outdir, f"{each_gene.replace('CDD:', '')}.faa"), 'w') as f1:
            SeqIO.write(unique_cdd_records, f1, format='fasta-2line')


def stats_cog(genome2genes, gene_ids):
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
    gene2genomes = {}
    for gene in gene_ids:
        _cache = [k for k, v in genome2genes.items() if v.get(gene, [])]
        # for genome, pdict in genome2genes.items():
        gene2genome_num[gene] = len(_cache)
        gene2genomes[gene] = _cache
    return gene_multi, gene_Ubiquity, gene2genome_num, gene2genomes


@click.command()
@click.option("-in_p", 'in_proteins', help='input directory which contains protein sequences file')
@click.option("-in_a", 'in_annotations', help="Actually output directory which contains annotations files during extraction")
@click.option("-s", "suffix", default='faa', help='suffix of protein files in `in_p`')
@click.option("-o", 'outdir', help="name of output directory")
@click.option("-evalue", 'evalue', default=1e-50, help="evalue for filtering out false-positive proteins. default is 1e-50 ")
@click.option("-gl", "genome_list", default=None,
              help="It will read 'selected_genomes.txt', please prepare the file, or indicate the alternative name or path. It could be None. If you provided, you could use it to subset the aln sequences by indicate names.")
@click.option("-ot", 'output_type', default='prot',help="prot(protein) or nucl(nucleotide)")
def main(in_proteins, suffix, in_annotations, outdir, evalue, genome_list,output_type):
    if genome_list is None:
        gids = []
    else:
        gids = open(genome_list).read().split('\n')
        gids = list(set([_ for _ in gids if _]))
    in_proteins = join(in_proteins, '*.' + suffix.strip('.'))
    protein_files = glob(in_proteins)
    if not protein_files:
        exit(f"error input proteins dir {in_proteins}")
    if not exists(in_annotations):
        os.makedirs(in_annotations)

    annotate_cog(protein_files, in_annotations)
    genome2cdd = parse_annotation(in_annotations, top_hit=True, evalue=evalue)
    if output_type.lower() in ['prot','protein']:
        write_cog(outdir, genome2cdd, in_proteins, genome_ids=gids, get_type='prot')
    elif output_type.lower() in ['nucl','nucleotide']:
        write_cog(outdir, genome2cdd, in_proteins, genome_ids=gids, get_type='nuc')
    else:
        raise IOError('wrong input of output_type')
    # write_cog(outdir + '_nuc', genome2cdd, in_proteins, genome_ids=gids, get_type='nuc')

    _subgenome2cdd = {k: v for k, v in genome2cdd.items() if k in set(gids)}
    gene_ids = set([_ for vl in genome2cdd.values() for _ in vl])
    gene_multi, gene_Ubiquity, gene2genome_num, gene2genomes = stats_cog(_subgenome2cdd, gene_ids)

    bb_g = [k for k, v in gene2genome_num.items() if v == len(gids)]
    if bb_g and gids:
        print(f"backbone genes is {str(bb_g)}")
    else:
        if genome_list:
            print("No backbone genes... all gene2genomes data could be reviewed at .. ")


if __name__ == "__main__":
    main()

# if __name__ == "__main__":
#     import sys

#     # usage :
#     # extract_cog.py 'raw_genome_proteins/*.faa' ./target_genes ./conserved_protein
#     if len(sys.argv) >= 2:
#         raw_proteins = sys.argv[1]
#         annotation_dir = sys.argv[2]
#         outdir = sys.argv[3]
#         gids = []
#     else:
#         raw_proteins = expanduser('~/data/nitrification_for/dating_for/raw_genome_proteins/*.faa')
#         annotation_dir = expanduser('~/data/nitrification_for/dating_for/target_genes_rpsblast')
#         outdir = expanduser('~/data/nitrification_for/dating_for/cog25_single')
#         gids = open(expanduser('~/data/nitrification_for/dating_for/bac120_annoate/remained_ids_fullv1.list')).read().split('\n')
