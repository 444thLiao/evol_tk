import sys
from os.path import *
import os
#sys.path.insert(0,dirname(dirname(dirname(dirname(__file__)))))
from subprocess import check_call
import click
from glob import glob
from tqdm import tqdm
import multiprocessing as mp
from Bio import SeqIO
command_template = 'mafft --maxiterate 1000 --genafpair --thread -1 {in_file} > {o_file} '
def run(args):
    unit_run(*args)
    
def unit_run(in_file,o_file):
    check_call(command_template.format(in_file=in_file,
                                       o_file=o_file), 
               shell=1)
    
def main(in_dir,odir,num_parellel,suffix='',new_suffix='',gids = None,force=False,**kwarg):
    suffix = suffix.strip('.')
    new_suffix = new_suffix.strip('.')
    if not exists(odir):
        os.makedirs(odir)
    if suffix:
        suffix = '.'+suffix
    file_list = glob(join(in_dir,f'*{suffix}'))
    if gids is not None:
        gids = set(gids)
        os.makedirs(join(odir,'tmp'),exist_ok=1)
        new_file_list = []
        for f in file_list:
            records = SeqIO.parse(f,format='fasta')
            records = [_ for _ in records if _.id in gids]
            n_f = join(odir,'tmp',basename(f))
            if not records:
                raise Exception('error not records')
            with open(n_f,'w') as f1:
                SeqIO.write(records,f1,format='fasta-2line')
            new_file_list.append(n_f)
        file_list = new_file_list[::]
    tqdm.write("start to process %s file with '%s' as suffix" % (len(file_list),suffix))
    params = []
    for in_file in tqdm(file_list):
        if new_suffix and suffix:
            ofile = join(odir,
                         basename(in_file).replace(suffix,
                                                   '.'+new_suffix))
        else:
            ofile = join(odir,
                         basename(in_file))
        if not exists(ofile) or force:
            params.append((in_file,ofile))
    with mp.Pool(processes=num_parellel) as tp:
        list(tp.imap(run,tqdm(params)))

@click.command()
@click.option('-i','indir')
@click.option('-o','odir')
@click.option('-s','suffix',default='faa')
@click.option('-ns','new_suffix',default='aln')
@click.option('-np','num_parellel',default=10)
@click.option("-gl", "genome_list", default=None, 
              help="it will read 'selected_genomes.txt', please prepare the file, or indicate the alternative name or path.")
@click.option('-f','force',help='overwrite?',default=False,required=False,is_flag=True)
def cli(indir,odir,num_parellel,suffix,new_suffix,genome_list,force):
    if genome_list is None:
        gids = None
    else:
        gids = open(genome_list).read().split('\n')
        gids = [_ for _ in gids if _]
    main(indir,odir,num_parellel,suffix,new_suffix,gids = gids,force=force)



if __name__ == "__main__":
    cli()