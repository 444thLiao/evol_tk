#!/usr/bin/env python

import os
import sys
import click
import time
from tqdm import tqdm
from os.path import realpath

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        print(dir)
        os.mkdir(dir)

zsh_path = os.popen('which zsh').read().strip('\n ')
job_directory = f"{os.getcwd()}/.job"
# Make top level directories
mkdir_p(job_directory)

def batch_iter(iter, batch_size):
    # generating batch according batch_size
    iter = list(iter)
    n_iter = []
    batch_d = 0
    for batch_u in range(0, len(iter), batch_size):
        if batch_u != 0:
            n_iter.append(iter[batch_d:batch_u])
        batch_d = batch_u
    n_iter.append(iter[batch_d: len(iter) + 1])
    return n_iter


def sbatch_all(cmds,reset_workdir=False,thread_per_tasks=1,fixed_cluster=''):
    count_ = 0
    for cmd in cmds:
        # cmds = open('./cmds').read().strip('\n').split('\n')
        # if reset_workdir:
        workdir = realpath(cmd.split(';')[0].strip().split(' ')[-1])
        cmd = cmd.split(';')[-1].strip()
        # cmd = cmd.replace()
        job_file = os.path.join(job_directory,f"job{count_}.job" )
        
        with open(job_file,'w') as fh:
            fh.writelines(f"#!{zsh_path}\n")
            fh.writelines(f"#SBATCH --job-name=job{count_}\n")
            fh.writelines(f"#SBATCH --cpus-per-task={thread_per_tasks}\n")
            fh.writelines(f"#SBATCH --output={job_directory}/job{count_}.out\n")
            if fixed_cluster:
                fh.writelines(f"#SBATCh -w {fixed_cluster} \n") 
            if reset_workdir:
                fh.writelines(f"#SBATCH --workdir={workdir}\n")
            fh.writelines(cmd)
        os.system("sbatch %s" % job_file)
        count_ += 1

@click.command()
@click.option("-i","infile",help="input file which stodge commands you want to batch")
@click.option("-w","reset_working",help="reseting the working directory. It mainly for mcmctree. ",
              default=False, required=False, is_flag=True,)
@click.option("-t","thread_per_tasks",default=None,help="Default would not set it. You could specify it to restrict the number of threads it used in each task.")
def cli(infile,reset_working,thread_per_tasks):
    cmds = open(infile).read().strip('\n').split('\n')
    sbatch_all(cmds,reset_workdir=reset_working,thread_per_tasks=int(thread_per_tasks))


if __name__ == "__main__":
    cli()
    
    from subprocess import check_call
    
    # cmd_file = sys.argv[1]
    # cmd_file = './cmds'
    # cmds = [_ for _ in open(cmd_file).read().strip('\n').split('\n')]

    # count_ = 0
    # for cmd in cmds:
    #     workdir = realpath(cmd.split(';')[0].strip().split(' ')[-1])
    #     cmd = cmd.split(';')[-1].strip()
    #     job_file = os.path.join(job_directory,f"mcmctree{count_}.job" )
        
    #     with open(job_file,'w') as fh:
    #         fh.writelines(f"#!{zsh_path}\n")
    #         fh.writelines(f"#SBATCH --job-name=mcmctree{count_}.job\n")
    #         fh.writelines(f"#SBATCH --cpus-per-task=1\n")
    #         fh.writelines(f"#SBATCH --workdir={workdir}\n")
    #         fh.writelines(cmd)
    #     os.system("sbatch %s" % job_file)
    #     # check_call(f"sbatch {job_file}",shell=1)
    #     count_ += 1

    # cmds = open('./cmds').read().split('\n')
    # from os.path import *
    # new_cmds = []
    # for cmd in tqdm(cmds):
    #     odir = cmd.strip().split(' ')[4]
    #     gbk_file = join(odir,basename(odir)+'.gbk')
    #     if not exists(gbk_file):
    #         new_cmds.append(cmd)

    # list_cmds = batch_iter(new_cmds,1500)

    # job_directory = './.job/'

    # for _cmds in list_cmds:
    #     count_ = 0
    #     for cmd in _cmds:
    #         # workdir = './'
    #         cmd = cmd.split(';')[-1]
    #         cmd = cmd.replace('--cpus 0','--cpus 10')
    #         job_file = os.path.join(job_directory, f"job_lth{count_}.job")

    #         with open(job_file, 'w') as fh:
    #             fh.writelines(f"#!{zsh_path}\n")
    #             fh.writelines(f"#SBATCH --job-name=job_lth{count_}.job\n")
    #             fh.writelines(f"#SBATCH --cpus-per-task=10\n")
    #             fh.writelines(f"#SBATCH --output={job_directory}/job_lth{count_}.out\n")
    #             #fh.writelines(f"#SBATCH --cluster-constraint=cl002\n")
                
    #             # fh.writelines(f"#SBATCH --error=.out/job_lth{count_}.err\n")
    #             # fh.writelines(f"#SBATCH --workdir={workdir}\n")

    #             fh.writelines(cmd)

    #         os.system("sbatch %s" % job_file)
    #         count_ += 1
    #     time.sleep(18000)  # 4hours later