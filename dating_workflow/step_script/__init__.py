from collections import defaultdict
from os.path import *
from subprocess import check_call

from Bio import SeqIO
from tqdm import tqdm
import os
from glob import glob
from api_tools.tk import run,refresh_tmp,get_files,get_genomes,get_tophit,parse_blastp,parse_hmmscan,convert_genome_ID,convert_genome_ID_rev

def process_path(path):
    if '~' in path:
        path = expanduser('path')
    if not '/' in path:
        path = './' + path
    if path.startswith('~'):
        path = expanduser(path)
    if path.startswith('.'):
        path = abspath(path)
    return path


def type_process(get_type):
    if get_type == 'nuc':
        suffix = 'ffn'
        final_suffix = 'ffn'
    elif get_type == 'prot':
        suffix = 'faa'
        final_suffix = 'faa'
    else:
        raise Exception('Unknown get_type %s' % str(get_type))
    return suffix, final_suffix


def try_get_file_from_formatted_dir(genome_name, 
                                    prokka_dir, 
                                    gene_file,
                                    suffix,
                                    gene2locus):
    genome2seq = {}
    collect_no_prokka_gids = []
    gfile = gene_file
    # try to find the prokka dir
    if prokka_dir is None:
        # Maybe it could reverse-seek the prokka directory.
        prokka_basedir = abspath(dirname(dirname(realpath(gfile))))
    else:
        # directly assign a directory
        prokka_basedir = abspath(realpath(prokka_dir))

    # the following format is for thliao scripts. (unchangeable)
    if suffix == 'faa':
        # get prot
        ori_file = gfile
    elif suffix == 'ffn':
        # get nucleotide of genes
        ori_file = f"{prokka_basedir}/{genome_name}/{genome_name}.{suffix}"
    else:
        raise IOError

    if exists(ori_file):
        _cache = {record.id: record
                  for record in SeqIO.parse(ori_file, format='fasta')}
        try:
            seq_set = {gene: [_cache[locus]
                       for (locus,evalue) in locus_list
                       if locus in _cache]
                   for gene, locus_list in gene2locus.items()}
        except:
            seq_set = {gene: [_cache[locus]
                       for locus in locus_list
                       if locus in _cache]
                   for gene, locus_list in gene2locus.items()}
        genome2seq[genome_name] = seq_set
    else:
        # not with prokka annotations
        print(genome_name, 'not annotated with prokka')
        collect_no_prokka_gids.append(genome_name)
        return genome2seq, collect_no_prokka_gids
    return genome2seq, collect_no_prokka_gids


# extract protein
def get_seq_and_write(outdir,
                      genome2cdd,
                      protein_files,
                      genome_ids=[],
                      get_type='prot',
                      _suffix='faa',
                      prokka_dir=None):
    """
    Extract the sequence and rename it with its genome name.
    It should be ensured that, for each gene/cdd, each genome contains only one locus.
    
    :param outdir: output directory
    :param genome2cdd: dict stored genome2cdd2locus
    :param protein_files: protein files 
    :param genome_ids: list of genome id
    :param _suffix: used to remove the suffix n protein files
    :param get_type: type of sequence want. nucl or prot
    :param prokka_dir: the directory of the output of prokka. Unless the protein file is a soft link.
    :return:
    """
    if not exists(outdir):
        os.makedirs(outdir)
    genome2seq = {}
    if not genome_ids:
        genome_ids = list(genome2cdd)
    gene_ids = set([_ for vl in genome2cdd.values() for _ in vl])
    suffix, final_suffix = type_process(get_type)

    tqdm.write('get sequence file')
    collect_no_prokka_gids = []
    for gene_file in tqdm(protein_files):
        genome_name = basename(gene_file).replace(f'.{_suffix}','')
        cdd2locus = genome2cdd[genome_name]
        if genome_name not in genome2cdd:
            tqdm.write(f'{genome_name} not in genome2cdd')
            continue
        _genome2seq, _collect_no_prokka_gids = try_get_file_from_formatted_dir(genome_name,
                                                                               prokka_dir,
                                                                               gene_file,
                                                                               suffix=suffix,
                                                                               gene2locus=cdd2locus)
        genome2seq.update(_genome2seq)
        collect_no_prokka_gids.extend(_collect_no_prokka_gids)
    if collect_no_prokka_gids:
        with open(join(outdir, 'no_prokka.gids'), 'w')  as f1:
            f1.write('\n'.join(collect_no_prokka_gids))
    # concat/output proteins
    tqdm.write('write out')
    for each_gene in tqdm(gene_ids):
        gene_records = []
        for gname, seq_dict in genome2seq.items():
            get_records = seq_dict.get(each_gene, [])
            for record in get_records:
                record.id = gname
            gene_records += get_records
        unique_cdd_records = []
        [unique_cdd_records.append(record)
         for record in gene_records
         if record.id not in [_.id for _ in unique_cdd_records]]

        with open(join(outdir, f"{each_gene.replace('CDD:', '')}.{final_suffix}"), 'w') as f1:
            SeqIO.write(unique_cdd_records, f1, format='fasta-2line')


def stats_cog(genome2genes, gene_ids):
    """
    get the basic stats of annotate genes. Such as:
    gene_multi: gene_name to the number of genomes containing multiple of it
    gene_Ubiquity: gene_name to the number of genomes containing it

    :param genome2genes:
    :param gene_ids:
    :return:
    """
    gene_multi = {g: 0 for g in gene_ids}
    for genome, pdict in genome2genes.items():
        for gene, seqs in pdict.items():
            if len(seqs) >= 2:
                gene_multi[gene] += 1
    gene_Ubiquity = {g: 0 for g in gene_ids}
    for genome, pdict in genome2genes.items():
        for gene, seqs in pdict.items():
            gene_Ubiquity[gene] += 1

    gene2genomes = {}
    for gene in gene_ids:
        _cache = [k for k, v in genome2genes.items() if v.get(gene, [])]
        # for genome, pdict in genome2genes.items():
        gene2genomes[gene] = _cache
    return gene_multi, gene_Ubiquity, gene2genomes


def write_out_stats(outdir, genome2genes, gene_ids):
    gene_multi, gene_Ubiquity, gene2genomes = stats_cog(genome2genes, gene_ids)
    if not exists(outdir):
        os.makedirs(outdir)
    with open(join(outdir, 'stats.tab'), 'w') as f1:
        f1.write("gene name\tnumber of genomes containing multiple of it\tnumber of genomes cotaining it\n")
        for g in gene_ids:
            f1.write(f"{g}\t{gene_multi[g]}\t{gene_Ubiquity[g]}\n")
    return gene_multi, gene_Ubiquity, gene2genomes

