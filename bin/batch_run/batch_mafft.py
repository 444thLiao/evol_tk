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

from dating_workflow.step_script import convert_genome_ID_rev

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


def main(in_dir, odir, num_parellel, suffix='', new_suffix='', gids=None, force=False, mode=default_mode, fix_refseq=False, removed_gene_list=None, not_add_prefix_ids=[],**kwarg):
    suffix = suffix.strip('.')
    new_suffix = new_suffix.strip('.')
    if not exists(odir):
        os.makedirs(odir)
    if suffix:
        suffix = '.' + suffix
    file_list = glob(join(in_dir, f'*{suffix}'))
    if gids is not None:
        gids = set(gids)
        os.makedirs(join(odir, 'tmp'), exist_ok=1)
        new_file_list = []
        tqdm.write('iterating files to collect with giving genome ids')
        for f in tqdm(file_list):
            records = SeqIO.parse(f, format='fasta')
            records = [_
                       for _ in records
                       if _.id in gids]
            if not records:
                records = SeqIO.parse(f, format='fasta')

                if not fix_refseq:
                    records = [_
                               for _ in records
                               if convert_genome_ID_rev(_.id.split('_')[0],not_add_prefix_ids=not_add_prefix_ids ) in gids]
                else:
                    gids = [_.split('_')[-1] for _ in gids]
                    records = [_
                               for _ in records
                               if convert_genome_ID_rev(_.id.split('_')[0], prefix='',not_add_prefix_ids=not_add_prefix_ids) in gids]
            n_f = join(odir, 'tmp', basename(f))
            if not records or len(records) == 1:
                print(f'failed records,for {f}, pass it')
                continue
            if removed_gene_list is not None:
                records = [_ for _ in records
                           if _.id not in removed_gene_list]
            with open(n_f, 'w') as f1:
                SeqIO.write(records, f1, format='fasta-2line')
            new_file_list.append(n_f)
        file_list = new_file_list[::]
    tqdm.write("start to process %s file with '%s' as suffix" % (len(file_list), suffix))
    params = []
    for in_file in tqdm(file_list):
        if new_suffix and suffix:
            ofile = join(odir,
                         basename(in_file).replace(suffix,
                                                   '.' + new_suffix))
        else:
            ofile = join(odir,
                         basename(in_file))
        if not exists(ofile) or force:
            params.append((in_file, ofile, mode))
    with mp.Pool(processes=num_parellel) as tp:
        r = list(tqdm(tp.imap(run, params), total=len(params)))


@click.command()
@click.option('-i', 'indir',help="input directory which stodge files with suffix. Normally, it should contain faa files directly.  ")
@click.option('-o', 'odir',help="output directory which stodge the ouput aln files.")
@click.option('-s', 'suffix', help='suffix for files', default='faa')
@click.option('-ns', 'new_suffix', default='aln')
@click.option('-np', 'num_parellel', default=10,help="number of parellel processes you want to run. ")
@click.option("-gl", "genome_list", default=None,
              help="It will read 'selected_genomes.txt', please prepare the file, or indicate the alternative name or path. / "
                   "It could be None. If you provided, you could use it to subset the aln sequences by indicate names.")
@click.option('-m', 'mode_mafft', default='ginsi',help="the mode of mafft you want to use. You could choose mafft, ginsi, einsi, linsi. You could find the detailed descriptions of them at the help of mafft. ")
@click.option('-f', 'force', help='overwrite?', default=False, required=False, is_flag=True)
@click.option('-fix_ref', 'fix_refseq', help='fix name of refseq?', default=False, required=False, is_flag=True)
@click.option('-rm_l', 'removed_gene_list', help='list of removed gene?')
@click.option('-not_add_prefix', 'not_add_prefix', help='provide a list of id which do not add prefix as others. ', default=None, required=False)
def cli(indir, odir, num_parellel, suffix, new_suffix, genome_list, force, mode_mafft, removed_gene_list, fix_refseq,not_add_prefix):
    if genome_list is None:
        gids = None
    else:
        gids = open(genome_list).read().split('\n')
        gids = set([_ for _ in gids if _])
    if removed_gene_list is not None:
        removed_gene_list = open(removed_gene_list).read().split('\n')
    if not_add_prefix is not None:
        not_add_prefix_ids = [_ for _ in open(not_add_prefix).read().split('\n') if _]
    else:
        not_add_prefix_ids = []
    main(indir, odir, num_parellel, suffix, new_suffix, gids=gids, force=force, mode=mode_mafft, removed_gene_list=removed_gene_list, fix_refseq=fix_refseq,not_add_prefix_ids=not_add_prefix_ids)


if __name__ == "__main__":
    cli()
