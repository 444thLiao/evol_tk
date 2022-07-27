"""
For easily perform mafft for some
"""

import multiprocessing as mp
import os
from glob import glob
from os.path import *
from subprocess import check_call

import click
from Bio import SeqIO
from tqdm import tqdm
from Bio.Seq import Seq
from dating_workflow.step_script import convert_genome_ID_rev,get_genomes,get_files,refresh_tmp

default_mode = 'ginsi'
command_template = '{mode} --anysymbol --thread -1 {in_file} > {o_file} '


def run(args):
    unit_run(*args)


def unit_run(in_file, o_file, mode):
    check_call(command_template.format(in_file=in_file,
                                       o_file=o_file,
                                       mode=mode),
               shell=True,
               stdout=open('/dev/null', 'w'),
               stderr=open('/dev/null', 'w'))

def remove_3rd(record,return_seq=False):
    # inplace remove the 3rd
    seq = str(record.seq)
    two_partitions = [seq[::3],seq[1::3]]
    final_seq = ''.join([''.join(_) for _ in zip(*two_partitions)])
    #final_seq = final_seq[:-2]
    if return_seq:
        return final_seq
    record.seq = Seq(final_seq)  # :-2 mean remove the lefted stop codon 
    return record

def main(in_dir, odir, num_parellel, suffix='', new_suffix='', 
         name2prefix=None, force=False, 
         mode=default_mode,
         gene_list=None,remove3rd=False):
    suffix = '.'+suffix.strip('.')
    new_suffix = '.'+new_suffix.strip('.')
    if not exists(odir):
        os.makedirs(odir)
    file_list = get_files(in_dir,suffix)
    if gene_list:
        file_list = [_ for _ in file_list if basename(_).replace(f"{suffix}",'') in gene_list]
        
    # prepare the tmp files which has been subsampled
    of2f_list = {}
    for f in file_list:
        n_f = join(odir, 'tmp', basename(f))
        if n_f in of2f_list:
            of2f_list[n_f].extend([f])
        else:
            of2f_list[n_f] = [f]
    if name2prefix is not None:
        all_prefix = [_ for k in name2prefix.values() for _ in k]
        refresh_tmp(join(odir, 'tmp'))
        new_file_list = []
        tqdm.write('iterating files to collect with giving genome ids')
        for n_f,f_list in tqdm(of2f_list.items()):
            # this part could be more slow to get the right selected records
            records = [_ for f in f_list for _ in SeqIO.parse(f, format='fasta')]
            records = [_
                       for _ in records
                       if _.id in all_prefix]
            if (not records) or (len(records) == 1):
                print(f'no available (or only one) record could be used in {f}')
                continue
            if remove3rd:
                records = [remove_3rd(r)
                           for r in records]
            with open(n_f, 'w') as f1:
                SeqIO.write(records, f1, format='fasta-2line')
            new_file_list.append(n_f)
    tqdm.write("start to process %s file with '%s' as suffix" % (len(new_file_list), suffix))
    params = []
    for in_file in tqdm(new_file_list):
        if new_suffix and suffix:
            ofile = join(odir,
                         basename(in_file).replace(suffix,
                                                   new_suffix))
        else:
            ofile = join(odir,
                         basename(in_file))
        if not exists(ofile) or force:
            params.append((in_file, ofile, mode))
    with mp.Pool(processes=num_parellel) as tp:
        r = list(tqdm(tp.imap(run, params), total=len(params)))


@click.command()
@click.option('-i', 'indir',help="input directory which stodge files with suffix. Normally, it should contain faa files directly. Multiple indir could be provided with comma-separated ")
@click.option('-o', 'odir',help="output directory which stodge the ouput aln files.")
@click.option('-s', 'suffix', help='suffix for files', default='faa')
@click.option('-ns', 'new_suffix', default='aln')
@click.option('-np', 'num_parellel', default=10,help="number of parellel processes you want to run. ")
@click.option("-gl", "genome_list", default=None,
              help="It will read 'selected_genomes.txt', please prepare the file, or indicate the alternative name or path. / "
                   "It could be None. If you provided, you could use it to subset the aln sequences by indicate names.")
@click.option('-m', 'mode_mafft', default='ginsi',help="the mode of mafft you want to use. You could choose mafft, ginsi, einsi, linsi. You could find the detailed descriptions of them at the help of mafft. ")
@click.option('-f', 'force', help='overwrite?', default=False, required=False, is_flag=True)
@click.option('-gene_l', 'gene_list', help='list of retained genes')
@click.option('-r3', 'remove3rd', help='remove 3rd codon', default=False, required=False, is_flag=True)
def cli(indir, odir, num_parellel, suffix, new_suffix, genome_list, force, mode_mafft, gene_list,remove3rd):
    
    name2prefix = get_genomes(genome_list,1)
    if gene_list is not None:
        gene_list = open(gene_list).read().split('\n')
    main(indir, odir, num_parellel, suffix, new_suffix, name2prefix=name2prefix, 
         force=force, mode=mode_mafft, gene_list=gene_list,remove3rd=remove3rd)


if __name__ == "__main__":
    cli()
