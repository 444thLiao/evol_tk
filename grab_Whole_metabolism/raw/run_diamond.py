from glob import glob
import click
import multiprocessing as mp
from subprocess import check_call
from os.path import join,dirname,basename,abspath

kegg_diamond_db = "/home-user/sswang/db/diamond/kegg/latest/kegg"
thread = 30

cmd_template = "diamond blastp -q {infile} -o {outfile} -d {kegg_diamond_db} -p {thread}"
def run_cmd(cmd):
    check_call(cmd,shell=True)

@click.command()
@click.option("-i","indir")
@click.option("-t","num_thread")
@click.option("-p","num_parallel")
def main(indir,num_thread,num_parallel):
    indir = abspath(indir)
    p_files = glob(join(indir,"*.faa"))
    if not p_files:
        raise Exception("No faa file in %s" % indir)
    
