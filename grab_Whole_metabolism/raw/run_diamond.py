from glob import glob
import click
import multiprocessing as mp
from subprocess import check_call
from os.path import join,dirname,basename,abspath
import os
from tqdm import tqdm

kegg_diamond_db = "/home-user/sswang/db/diamond/kegg/latest/kegg"
thread = 30

cmd_template = "diamond blastp --quiet -q {infile} -o {outfile} -d {kegg_diamond_db} -p {thread}"
def run_cmd(cmd):
    check_call(cmd,shell=True)

# @click.command()
# @click.option("-i","indir")
# @click.option("-o","outdir")
# @click.option("-t","num_thread",default=30)
# @click.option("-p","num_parallel",default=2)
def main(indir,outdir,num_thread,num_parallel):
    indir = abspath(indir)
    p_files = glob(join(indir,'**',"*.faa"))
    if not p_files:
        raise Exception("No faa file in %s" % indir)
    if not os.path.exists(abspath(outdir)):
        os.makedirs(outdir,exist_ok=True)

    cmd_list = []
    for in_f in p_files:
        name = in_f.rpartition('.')[0]
        cmd = cmd_template.format(infile=in_f,
                                  outfile=join(outdir,name+'.diamond_out'),
                                  kegg_diamond_db=kegg_diamond_db,
                                  thread=num_thread)
        cmd_list.append(cmd)
    with mp.Pool(processes=num_parallel) as tp:
        for r in tqdm(tp.imap(run_cmd,cmd_list),total=len(cmd_list)):
            pass

if __name__ == '__main__':
    main()
