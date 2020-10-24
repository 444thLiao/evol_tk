import multiprocessing as mp
import os
import pickle
from collections import defaultdict
from glob import glob
from os.path import *
import hashlib
import click
from tqdm import tqdm

from dating_workflow.step_script import parse_hmmscan, run, get_seq_and_write, write_out_stats

HOME = os.getenv("HOME")
pfam_db = f'{HOME}/data/protein_db/bac120/Pfam.v32.sub6.hmm'
tigfam_db = f'{HOME}/data/protein_db/bac120/TIGRFAMv14_sub114.hmm'

__file__ = f'{HOME}/script/evolution_relative/dating_workflow/step_script/extrat_bac120.py'
bac120_list = join(dirname(__file__), 'bac120.tsv')
id_list = [row.split('\t')[0] for row in open(bac120_list) if row]
id_list = id_list[1:]  # remove first row
pfam_ids = [_.strip() for _ in id_list if _.startswith('PF0')]
tigfam_ids = [_.strip() for _ in id_list if _.startswith('TIGR')]


# ABOVE is the default setting for luolab server.

def annotate_bac120(protein_files, odir, db_id='pfam',cpu=10,num_p=5):
    params = []
    if not exists(odir):
        os.makedirs(odir)
    
    size0_pfiles = []
    hmmscan = '`which hmmscan`'
    for pfile in protein_files:
        gname = basename(pfile).replace('.faa', '')
        if getsize(pfile) !=0:
            size0_pfiles.append(pfile)
            continue
        if db_id == 'pfam':
            ofile = f'{odir}/PFAM/{gname}.out'
            cmd = f"{hmmscan} --tblout {ofile} --acc --noali --notextw --cpu {cpu} {pfam_db} {pfile}"
        elif db_id == 'tigrfam':
            ofile = f'{odir}/TIGRFAM/{gname}.out'
            cmd = f"{hmmscan} --tblout {ofile} --acc --noali --notextw --cpu {cpu} {tigfam_db} {pfile}"
        else:
            raise SyntaxError('unknown %s' % db_id)
        if not exists(ofile):
            if not exists(dirname(ofile)):
                os.makedirs(dirname(ofile))
            params.append(cmd)
            # check_call(cmd, shell=1)
    # print(params)
    with mp.Pool(processes=num_p) as tp:
        r = list(tqdm(tp.imap(run, params), total=len(params)))


def parse_annotation(odir, top_hit=False, evalue=1e-50):
    # for cdd
    _cdd_match_ids = pfam_ids
    genome2annotate = defaultdict(lambda: defaultdict(list))

    # cdd annotations
    
    cdd_anno_files = glob(join(odir, 'PFAM', '*.out'))
    # tigrfam annotations
    tigrfam_anno_files = glob(join(odir, 'TIGRFAM', '*.out'))
    # add cache to avoid iterate it again and again
    # m = hashlib.md5()
    # hash_str = m.update(tuple(sorted(tigrfam_anno_files + cdd_anno_files + [str(top_hit)])))
    # cache_file = join(odir, f'.tmp{hash_str}')
    # if exists(cache_file):
    #     genome2annotate = pickle.load(open(cache_file, 'rb'))
    #     genome2annotate = dict(genome2annotate)
    #     return genome2annotate
    tqdm.write('start to read/parse output files (cdd and tigrfam)')
    for ofile in tqdm(tigrfam_anno_files + cdd_anno_files):
        gname = basename(ofile).replace('.out', '')
        locus_dict = parse_hmmscan(ofile=ofile,
                                   top_hit=top_hit,
                                   filter_evalue=evalue,
                                   gene_pos=1)
        genome2annotate[gname].update(locus_dict)
        
    genome2annotate = dict(genome2annotate)
    # if not exists(cache_file):
    #     with open(cache_file, 'wb') as f1:
    #         pickle.dump(genome2annotate, f1)
    return genome2annotate


@click.command()
@click.option("-in_p", 'in_proteins', help='input directory which contains protein sequences file')
@click.option("-in_a", 'in_annotations', help="Actually output directory which contains annotations files during extraction")
@click.option("-s", "suffix", default='faa', help='suffix of protein files in `in_p`')
@click.option("-o", 'outdir', help="name of output directory",default=None,required=False)
@click.option("-evalue", 'evalue', default=1e-50, help="evalue for filtering out false-positive proteins. default is 1e-50 ")
@click.option("-ot", 'output_type', default='prot', help="prot(protein) or nucl(nucleotide)")
@click.option("-gl", "genome_list", default=None,
              help="It will read 'selected_genomes.txt', please prepare the file, or indicate the alternative name or path. It could be None. If you provided, you could use it to subset the aln sequences by indicate names.")
@click.option("-pd", 'prokka_dir', default=None,
              help="directory which restore the output of each genome. acceptable directory should contain  prokka_dir/genome_id/genome_id.faa and ffn (for nucleotide). ")
@click.option("-p_a", 'pass_annotation', is_flag=True,default=False,
              help="skip the annotation parts")
@click.option("-a", 'annotation_only', is_flag=True,default=False,
              help="only run the annotation parts")
def main(in_proteins, suffix, in_annotations, outdir, evalue, genome_list, output_type, prokka_dir,pass_annotation,annotation_only):
    if genome_list is None:
        gids = []
    else:
        gids = open(genome_list).read().split('\n')
        gids = list(set([_ for _ in gids if _]))
    protein_files = glob(join(in_proteins, '*.' + suffix.strip('.')))
    # gids = []
    if not protein_files:
        exit(f"error input proteins dir {in_proteins}")
    tqdm.write("Annotating these proteins, it only run once.. For tigrfam and pfam.")
    
    if not pass_annotation:
        annotate_bac120(protein_files, in_annotations, db_id='tigrfam')
        annotate_bac120(protein_files, in_annotations, db_id='pfam')
    if annotation_only:
        exit(f"finish annotation")
        
    tqdm.write("Parsing the annotation results...")
    genome2genes = parse_annotation(in_annotations, top_hit=False)
    gene_ids = pfam_ids + tigfam_ids

    _subgenome2cdd = {k: v for k, v in genome2genes.items() if k in set(gids)}
    write_out_stats(outdir,_subgenome2cdd, gene_ids)

    genome2genes = parse_annotation(
        in_annotations, top_hit=True, evalue=evalue)
    if gids:
        _subgenome2cdd = {k: v for k, v in genome2genes.items() if k in set(gids)}
    else:
        _subgenome2cdd = genome2genes.copy()
    if output_type.lower() in ['prot', 'protein']:
        get_seq_and_write(outdir, _subgenome2cdd, in_proteins, genome_ids=gids, get_type='prot', prokka_dir=prokka_dir)
    elif output_type.lower() in ['nucl', 'nucleotide']:
        get_seq_and_write(outdir, _subgenome2cdd, in_proteins, genome_ids=gids, get_type='nuc', prokka_dir=prokka_dir)
    else:
        raise IOError('wrong input of output_type')


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
