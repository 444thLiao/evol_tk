from collections import defaultdict
from os.path import *
from subprocess import check_call

from Bio import SeqIO
from tqdm import tqdm
import os

def run(cmd):
    check_call(cmd,
               shell=True,
               stdout=open('/dev/null', 'w'),
               )


def get_tophit(gid2locus, top_hit):
    if top_hit:
        gid2locus = {k: sorted(v,
                               key=lambda x: x[1])
                     for k, v in gid2locus.items()}
        gid2locus = {k: [v[0][0]]
        if v else []
                     for k, v in gid2locus.items()}
    else:
        gid2locus = {k: [_[0] for _ in v]
        if v else []
                     for k, v in gid2locus.items()}
    return gid2locus


def parse_blastp(ofile, match_ids=[], filter_evalue=1e-3, top_hit=False):
    if not match_ids:
        gid2locus = defaultdict(list)
    else:
        gid2locus = {k: [] for k in match_ids}
    for row in open(ofile, 'r'):
        sep_v = row.split('\t')
        locus = sep_v[0]
        evalue = float(sep_v[10])
        if filter_evalue and evalue > filter_evalue:
            continue
        if sep_v[1] in match_ids:
            gid2locus[sep_v[1]].append((locus, evalue))
        if not match_ids:
            gid2locus[sep_v[1]].append((locus, evalue))
    gid2locus = get_tophit(gid2locus, top_hit=top_hit)
    return gid2locus


def parse_hmmscan(ofile, filter_evalue=1e-20, top_hit=False, gene_pos=0):
    gid2locus = defaultdict(list)

    for row in open(ofile, 'r'):
        if row.startswith('#'):
            continue
        r = row.split(' ')
        r = [_ for _ in r if _]

        gene_id = r[gene_pos]
        locus_tag = r[2]
        evalue = float(r[4])
        if filter_evalue and evalue > filter_evalue:
            continue
        gid2locus[gene_id].append((locus_tag, evalue))
    if top_hit:
        gid2locus = get_tophit(gid2locus, top_hit=top_hit)
    return gid2locus


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


def try_get_file_from_formatted_dir(genome_name, prokka_dir, suffix, protein_dir,
                                    gene2locus):
    genome2seq = {}
    collect_no_prokka_gids = []
    gfile = f'{protein_dir}/{genome_name}.faa'
    # try to find the prokka dir
    if prokka_dir is None:
        # Maybe it could reverse-seek the prokka directory.
        prokka_basedir = abspath(dirname(dirname(realpath(gfile))))
    else:
        # directly assign a directory
        prokka_basedir = abspath(realpath(prokka_dir))

    # the following format is came from thliao scripts. (unchangeable)
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
                      protein_dir,
                      genome_ids=[],
                      get_type='prot',
                      prokka_dir=None):
    """

    :param outdir: output directory
    :param genome2cdd: dict stored genome2cdd2locus
    :param used_protein_file: the protein file used in annotation
    :param genome_ids: genome id
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
    for genome_name in tqdm(genome_ids):
        cdd2locus = genome2cdd[genome_name]

        _genome2seq, _collect_no_prokka_gids = try_get_file_from_formatted_dir(genome_name,
                                                                               prokka_dir,
                                                                               suffix,
                                                                               protein_dir,
                                                                               cdd2locus)
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
                record.name = gname
            gene_records += get_records
        unique_cdd_records = []
        [unique_cdd_records.append(record)
         for record in gene_records
         if record.id not in [_.id
                              for _ in unique_cdd_records]]

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


# two function for dating workflow(formatting the name of assembly genomes)
def convert_genome_ID(genome_ID):
    # for GCA_900078535.2
    # it will return 900078535v2
    if isinstance(genome_ID, str) and genome_ID.startswith('GC'):
        return genome_ID.split('_')[-1].replace('.', 'v')
    else:
        return genome_ID


def convert_genome_ID_rev(locus_ID, prefix='GCA_',not_add_prefix_ids=[]):
    # for 900078535v2
    # it will return prefix + 900078535.2
    if locus_ID in not_add_prefix_ids:
        return locus_ID
    if '|' in str(locus_ID):
        # other labmate used
        genome_name = locus_ID.partition('|')[0]
        return genome_name

    if isinstance(locus_ID, str) and not locus_ID.startswith('GC'):
        if '_' in locus_ID:
            # tianhua version, it won't contain |
            locus_ID = locus_ID.partition('_')[0]
            if locus_ID in not_add_prefix_ids:
                return locus_ID
            else:
                return prefix + locus_ID.replace('v', '.')
        else:
            return prefix + locus_ID.replace('v', '.')
    else:
        return locus_ID
