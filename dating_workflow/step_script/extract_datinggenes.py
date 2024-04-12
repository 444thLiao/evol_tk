"""
Input bacterial genomes protein files


Annotate them and generate  candidate files for alignment, trimal, and concatenated (including eukaryotic genomess)

"""
import copy
import hashlib
import multiprocessing as mp
import os
import pickle
from collections import *
from collections import defaultdict
from glob import glob
from os.path import *

import click
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from dating_workflow.step_script import (get_files, get_genomes,
                                         get_seq_and_write, parse_blastp,
                                         parse_hmmscan, run, write_out_stats)

####################################################################
#### COG25
####################


#### nucleus
####################


#### MITO
####################
dataset2usage = {"Oryza sativa": "both",
                 "Arabidopsis thaliana": "both",
                 "Physcomitrella patens": "both",
                 "Ostreococcus tauri": "both",
                 "Chondrus crispus": "both",
                 "Porphyra umbilicalis": "both",
                 "Cyanidioschyzon merolae": "both",
                 "Cyanophora paradoxa": "both",
                 "Phytophthora infestans": "both",
                 "Andalucia godoyi": "both",
                 "Jakoba libera": "mito",
                 "Jakoba bahamiensis": "mito",
                 "Reclinomonas americana": "mito",
                 "Histiona aroides": "mito",
                 "Seculamonas ecuadoriensis": "mito",
                 "Malawimonas jakobiformis": "mito",
                 "Thalassiosira pseudonana": "nuclear",
                 "Symbiodinium minutum": "nuclear",
                 "Paramecium tetraurelia": "nuclear",
                 "Oxytricha trifallax": "nuclear",
                 "Reticulomyxa filosa": "nuclear",
                 "Elphidium margaritaceum": "nuclear",
                 "Homo sapiens": "nuclear",
                 "Gallus gallus": "nuclear",
                 "Branchiostoma floridae": "nuclear",
                 "Amphimedon queenslandica": "nuclear",
                 "Candida albicans": "nuclear",
                 "Ustilago maydis": "nuclear",
                 "Pleurotus ostreatus": "nuclear",
                 "Spizellomyces punctatus": "nuclear",
                 "Dictyostelium discoideum": "nuclear",
                 "Polysphondylium pallidum": "nuclear",
                 "Acanthamoeba castellanii": "nuclear", }

mito24_proteins_reference = "/mnt/ivy/thliao/db/dating_for/mito_faa"
mito24_nucleotide_reference = "/mnt/ivy/thliao/db/dating_for/mito_ffn"
# copy files from  "/home-user/sswang/project/Mito/results/euk_tree/date/acc/plant-nobc/pep-sing/"
nuc22_nucleotide_reference = "/home-user/sswang/project/Mito/results/euk_tree/date/acc-cds/M9P-nobc/pep-sing/"
# full length of genes instead of 3rd trimmed sequences
nuc22_protein_reference ="/home-user/sswang/project/Mito/results/euk_tree/date/acc/full.M9P/pep-sing/"

# Noted that
# these should be only taken as reference for annotations
# if you want to get the mito24 genes of eukaryotes, you should see '/home-user/sswang/project/Mito/data/Pisani-cds_retrieval/euk-mito/cds-query/'. Under this folder, each organism has its own mito24 protein sequences.
# For nuc22, you should see /home-user/sswang/project/Mito/data/cds_retrieval/all/cds-query/


#### plastid
####################
plastid97_prot_Blastdb = f"/mnt/ivy/thliao/db/dating_for/plastid97_prot/plastid"
plastid97_prot_faa = f"/mnt/ivy/thliao/db/dating_for/plastid97_prot/plastid.faa"

EGT_prot_faa = f"/mnt/ivy/thliao/db/dating_for/plastidEGT_prot/EGT_nucleus.faa"
EGT_prot_Blastdb = f"/mnt/ivy/thliao/db/dating_for/plastidEGT_prot/EGT_nucleus"

####################################################################
######## Function
list_db = ["Plastid","Plastid_nucleus","mito24","nuc22","cog25","gomez19",
           "all"]
db2blastdb = {'plastid':plastid97_prot_Blastdb,
              'Plastid_nucleus':EGT_prot_Blastdb,
              }
db2hmmdb = {}


def anno_blastp(faa_list, _outdir,db,num_p=5,redo=False):
    cmds = []
    size0_pfiles = []
    if not exists(_outdir):
        os.makedirs(_outdir)
    for in_p in tqdm(faa_list):
        if getsize(in_p) ==0:
            size0_pfiles.append(in_p)
            continue
        aid = basename(in_p).rpartition('.')[0]
        ofile = f"{_outdir}/{aid}.tab"
        cmd = f"blastp -query '{in_p}' -db {db} -outfmt 6 -out '{ofile}' -evalue 1e-10 -max_hsps 20 -num_threads 20"
        if not exists(ofile) or getsize(ofile) == 0 or redo:
            cmds.append(cmd)
    tqdm.write(f'{len(size0_pfiles)} files are empty.')
    with mp.Pool(processes=num_p) as tp:
        list(tqdm(tp.imap(run, cmds), total=len(cmds)))

def parse_annotation(anno_files, odir=None,evalue=1e-10,mode='hmm',
                     db2length={},genome_ids = []):
    if not exists(odir) and odir is not None:
        os.system(f"mkdir -p {odir}")
    # else:
    #     os.system(f"rm {odir}/*.faa")
    genome2annotate = defaultdict(dict)
    for ofile in tqdm(anno_files):
        gname = basename(ofile).replace('.tab', '')
        if mode == 'hmm':
            locus2gene,gene2locus = parse_hmmscan(ofile=ofile,
                                   top_hit=False,
                                   filter_evalue=evalue,
                                   _pos=[2,0])

        elif mode == 'blastp':
            locus2gene,gene2locus = parse_blastp(ofile,
                                      filter_evalue=evalue,
                                      top_hit=False,
                                      seq2length=db2length,
                                      pos_location=[8,9,1,0,1],
                                       cov=60)
        # gene2locus_l = defaultdict(list)
        # for l,gene in locus_dict.items():
        #     gene2locus_l[gene[0]].append(l) # only mapped to a single gene
        # locus_dict = {k:list(v) for k,v in gene2locus_l.items()}
        gene2locus = {k:[_[0] for _ in v] for k,v in gene2locus.items()}
        genome2annotate[gname].update(gene2locus)

    df = pd.DataFrame.from_dict(genome2annotate)
    bin_df = df.applymap(lambda x:len(x) if str(x)!='nan' else 0)
    alllocus_df = df.applymap(lambda x:';'.join(x) if str(x)!='nan' else 'NA')
    finallocus_df = df.applymap(lambda x:x[0] if str(x)!='nan' else 'NA')
    if genome_ids:
        bin_df = bin_df.reindex(columns=genome_ids)
        alllocus_df = alllocus_df.reindex(columns=genome_ids)
        finallocus_df = finallocus_df.reindex(columns=genome_ids)
    if bin_df.shape[0] == 0:
        print(f"something wrong. No genomes remained. ")
    stats_df = get_stats(bin_df)
    
    if odir is not None:
        alllocus_df.to_csv(f'{odir}/matched_locus.tsv',sep='\t',index=1,index_label='Genes')
        finallocus_df.to_csv(f'{odir}/Tophit_locus.tsv',sep='\t',index=1,index_label='Genes')
        stats_df.to_csv(f'{odir}/stats.tsv',sep='\t',index=1,index_label='Genes')
    return alllocus_df,stats_df,finallocus_df

def get_stats(bin_df):
    """
    rows are genes. columns are genomes
    the vals should be the concatenated string of all matched locuse separated by ';'

    Args:
        df (_type_): _description_

    Returns:
        _type_: _description_
    """
    s_df = pd.DataFrame(index=bin_df.index,columns=['# Genomes (total)','# Genes','Mean copy number','Presence ratio'])
    s_df.loc[:,'# Genomes (total)'] = len(bin_df.columns)
    s_df.loc[:,'# Genes'] = bin_df.sum(1)
    s_df.loc[:,'Mean copy number'] = bin_df.sum(1)/(bin_df!=0).sum(1)
    s_df.loc[:,'Presence ratio'] = (bin_df!=0).sum(1)/bin_df.shape[1]
    return s_df


def extract_faa(alllocus_df,odir,genome2locus2seq):
    if not exists(odir):
        os.makedirs(odir)
    for gene,row in alllocus_df.iterrows():
        seqs = []
        for genome,locus_str in row.to_dict().items():
            if locus_str in ['','NA']:
                continue
            for l in locus_str.split(';'):
                seq = copy.deepcopy(genome2locus2seq[genome][l])
                seqs.append(seq)
        with open(f"{odir}/{gene}.faa",'w') as f1:
            SeqIO.write(seqs,f1,'fasta-2line')


####################################################################
######## ! preparing necessary data

def get_db2length(db):
    db2length = {}
    if db == 'Plastid':
        for r in SeqIO.parse(plastid97_prot_faa,'fasta'):
            db2length[r.id] = len(r.seq)
    elif db == 'Plastid_nucleus':
        for r in SeqIO.parse(f"/mnt/ivy/thliao/db/dating_for/plastidEGT_prot/EGT_nucleus.faa",'fasta'):
            db2length[r.id] = len(r.seq)
    elif db == 'mito24':
        for r in SeqIO.parse('/mnt/home-backup/thliao/cyano/analysis/nfixer/geneset/EGT_nucleus.faa','fasta'):
            db2length[r.id] = len(r.seq)
    return db2length




@click.command()
@click.option("-in_p", 'in_proteins', help='input directory which contains protein sequences file')
@click.option("-in_a", 'in_annotations', help="Actually output directory which contains annotations files during extraction")
@click.option("-s", "suffix", default='faa', help='suffix of protein files in `in_p`')
@click.option('-db','--db',
              type=click.Choice(list_db, case_sensitive=False))
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
         output_type,db, prokka_dir,pass_annotation,annotation_only,num_parellel):
    gids = get_genomes(genome_list,True)
    protein_files = get_files(in_proteins,suffix)
    if gids:
        protein_files = [_ for _ in protein_files if basename(_).replace(f'.{suffix}','') in gids]
    if not protein_files:
        exit(f"Wrong input proteins dir '{in_proteins}' since no protein files were found")
    if not exists(in_annotations):
        os.makedirs(in_annotations)
    in_annotations = realpath(in_annotations)
    if not pass_annotation:
        if db in db2blastdb:
            dbpath = db2blastdb[db]
            search_type = 'blastp'
            anno_blastp(protein_files, in_annotations,
                    db=dbpath,num_p=num_parellel,)
        
        
    if annotation_only:
        exit(f"finish annotation")
    db2length = get_db2length(db)
    anno_files = glob(in_annotations+'/*.tab')
    alllocus_df, stats_df, finallocus_df = parse_annotation(anno_files,
                                                          odir=in_annotations,db2length=db2length,mode=search_type,genome_ids=gids)
    sub_statsdf = stats_df.loc[(stats_df.loc[:,'Mean copy number']<=2) & (stats_df.loc[:,'Presence ratio']>=0.3),:]
    print(f"Annotatable genes: {stats_df.shape[0]}; Remaining genes: {sub_statsdf.shape[0]}",)
    tqdm.write('read original protein files and then extract')
    genome2locus2seq = {basename(faa).replace(f'.{suffix}',''):{r.id:r 
                                                                for r in SeqIO.parse(faa,'fasta')}
                    for faa in tqdm(protein_files)
                    if faa.split('/')[-1].replace(f'.{suffix}','') in gids}
    extract_faa(finallocus_df.loc[sub_statsdf.index,list(gids) ] ,
                outdir,
                genome2locus2seq)


if __name__ == '__main__':

    main()


    # cd /mnt/storage3/thliao/project/AOB/extend_MOB/
    # python ~/script/evol_tk/dating_workflow/step_script/extract_cog25.py -i ./raw_proteins/*.faa -in_a ./annotations/Plastid -db plastid 




