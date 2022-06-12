"""
extract 25 proteins for dating analysis
"""
import multiprocessing as mp
import os
import pickle
from collections import defaultdict
from glob import glob
from os.path import *
import hashlib
import click
from tqdm import tqdm

from dating_workflow.step_script import parse_blastp, run, get_seq_and_write, write_out_stats,get_genomes,get_files

HOME = os.getenv("HOME")
resource_dir = f"{HOME}/data/protein_db/dating_resource"
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

def annotate_cog(protein_file_list, cog_out_dir,num_p=5,suffix='.faa'):
    params = []
    size0_pfiles = []
    for f in protein_file_list:
        gname = basename(f).replace(suffix, '')
        # for cdd
        if getsize(f) ==0:
            size0_pfiles.append(f)
            continue
        ofile = f'{cog_out_dir}/{gname}.out'
        cmd = f"`which rpsblast` -query {f} -db {cog_db} -max_target_seqs 10 -num_threads 10 -outfmt 6 -evalue 1e-3  -out {ofile}"
        if (not exists(ofile)) or (getsize(ofile)==0):
            # in case manually interrupted
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
    tqdm.write(f'{len(size0_pfiles)} files are empty.')
    with mp.Pool(processes=num_p) as tp:
        list(tqdm(tp.imap(run, params), total=len(params)))


def parse_annotation(cog_out_dir, top_hit=False, evalue=1e-3):
    cog_out_dir = realpath(cog_out_dir) # otherwise the following hashlib might be differnt
    # for cdd
    # _cdd_match_ids = set([_ for vl in cdd_num.values() for _ in vl])
    genome2cdd = defaultdict(lambda: defaultdict(list))

    # cdd annotations
    tqdm.write('start to read/parse output files')
    cdd_anno_files = glob(join(cog_out_dir, '*.out'))

    t = ''.join(sorted(cdd_anno_files + [str(top_hit)] + [str(evalue)] ))
    m = hashlib.md5(t.encode())
    hash_str = m.hexdigest()
    cache_file = join(cog_out_dir, f'.tmp{hash_str}')
    if exists(cache_file):
        genome2cdd = pickle.load(open(cache_file, 'rb'))
        genome2cdd = dict(genome2cdd)
        return genome2cdd

    for ofile in tqdm(cdd_anno_files):
        gname = basename(ofile).replace('.out', '')
        locus2gene,gene2list_locus = parse_blastp(ofile=ofile,
                                                  match_ids=[], 
                                                  pos_location = [8,9,0,0,1],
                                                  top_hit=top_hit,  
                                                  filter_evalue=evalue)
        genome2cdd[gname].update({k:[v] for k,v in gene2list_locus.items()})
    genome2cdd = dict(genome2cdd)  
    
    # change it into normal dict in order to pickle it
    # if not exists(cache_file):
    os.system(f"find {dirname(cache_file)} -mtime +2 -name '.tmp*' -delete ")  # delete 2days ago cache
    with open(cache_file, 'wb') as f1:
        pickle.dump(genome2cdd, f1)
    return genome2cdd


@click.command()
@click.option("-in_p", 'in_proteins', help='input directory which contains protein sequences file')
@click.option("-in_a", 'in_annotations', help="Actually output directory which contains annotations files during extraction")
@click.option("-s", "suffix", default='faa', help='suffix of protein files in `in_p`')
@click.option("-o", 'outdir', help="name of output directory",default=None,required=False)
@click.option("-evalue", 'evalue', default=1e-20, help="evalue for filtering out false-positive proteins. default is 1e-20 ")
@click.option("-gl", "genome_list", default=None,
              help="It will read 'selected_genomes.txt', please prepare the file, or indicate the alternative name or path. It could be None. If you provided, you could use it to subset the aln sequences by indicate names.")
@click.option("-ot", 'output_type', default='prot', help="prot(protein) or nucl(nucleotide)")
@click.option("-pd", 'prokka_dir', default=None,
              help="directory which restore the output of each genome. acceptable directory should contain  prokka_dir/genome_id/genome_id.faa and ffn (for nucleotide). ")
@click.option("-p_a", 'pass_annotation', is_flag=True,default=False,
              help="skip the annotation parts")
@click.option("-a", 'annotation_only', is_flag=True,default=False,
              help="only run the annotation parts")
@click.option('-np', 'num_parellel', default=5, help="num of processes could be parellel.. default is 10")
def main(in_proteins, suffix, in_annotations, outdir, evalue, genome_list, 
         output_type, prokka_dir,pass_annotation,annotation_only,num_parellel):
    # if genome_list is None:
    #     gids = []
    # else:
    #     gids = open(genome_list).read().split('\n')
    #     gids = list(set([_ for _ in gids if _]))
    gids = get_genomes(genome_list,True)
    protein_files = get_files(in_proteins,suffix)
    if gids:
        protein_files = [_ for _ in protein_files if basename(_).replace(f'.{suffix}','') in gids]
    if not protein_files:
        exit(f"error input proteins dir {in_proteins} since no wanted files were found")
    if not exists(in_annotations):
        os.makedirs(in_annotations)
    if not pass_annotation:
        annotate_cog(protein_files, in_annotations,num_p=num_parellel,suffix=f'.{suffix}')
    if annotation_only:
        exit(f"finish annotation")
    genome2cdd = parse_annotation(in_annotations, top_hit=True, evalue=evalue)
    if gids:
        _subgenome2cdd = {k: v for k, v in genome2cdd.items() if k in set(gids)}
    else:
        _subgenome2cdd = genome2cdd.copy()
    if output_type.lower() in ['prot', 'protein']:
        get_seq_and_write(outdir, 
                          genome2cdd, 
                          protein_files, 
                           _suffix=suffix,
                          get_type='prot', 
                          prokka_dir=prokka_dir)
        # output sequence would use genome name as its sequence id
    elif output_type.lower() in ['nucl', 'nucleotide']:
        get_seq_and_write(outdir, genome2cdd, protein_files,_suffix=suffix, get_type='nuc', prokka_dir=prokka_dir)
    else:
        raise IOError('wrong input of output_type')

    _subgenome2cdd = {k: v for k, v in genome2cdd.items() if k in set(gids)}
    gene_ids = set([_ for vl in genome2cdd.values() for _ in vl])
    gene_multi, gene_Ubiquity, gene2genomes = write_out_stats(outdir,_subgenome2cdd, gene_ids)

    bb_g = [k
            for k, v in gene_Ubiquity.items()
            if v == len(gids)]
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
