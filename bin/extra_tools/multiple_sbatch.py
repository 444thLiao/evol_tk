#!/usr/bin/env python

import os
import sys
import click
import time
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

@click.command()
@click.option("-i","infile",help="input file which stodge commands you want to batch")
@click.option("")
def cli():
    pass


if __name__ == "__main__":

    cmd_files = sys.argv[1]
    cmds=[_ for _ in open(cmd_files).read().split('\n')]

    count_ = 0
    for cmd in cmds:
        workdir = cmd.split(';')[0].strip().split(' ')[-1]
        cmd = cmd.split(';')[-1]
        job_file = os.path.join(job_directory,f"job_lth{count_}.job" )


        with open(job_file,'w') as fh:
            fh.writelines(f"#!{zsh_path}\n")
            fh.writelines(f"#SBATCH --job-name=job_lth{count_}.job\n")
            # fh.writelines(f"#SBATCH --output=.out/job_lth{count_}.out\n")
            # fh.writelines(f"#SBATCH --error=.out/job_lth{count_}.err\n")
            fh.writelines(f"#SBATCH --workdir={workdir}\n")
            
            fh.writelines(cmd)

        os.system("sbatch %s" %job_file)
        count_ += 1

    cmds = open('./cmds').read().split('\n')
    from os.path import *
    new_cmds = []
    for cmd in cmds:
        odir = cmd.strip().split(' ')[4]
        gbk_file = join(odir,basename(odir)+'.gbk')
        if not exists(gbk_file):
            new_cmds.append(cmd)

    list_cmds = batch_iter(new_cmds,1500)

    job_directory = './.job/'

    for _cmds in list_cmds:
        count_ = 0
        for cmd in _cmds:
            # workdir = './'
            cmd = cmd.split(';')[-1]
            # cmd = cmd.replace('--cpus 30','--cpus 30')
            job_file = os.path.join(job_directory, f"job_lth{count_}.job")

            with open(job_file, 'w') as fh:
                fh.writelines(f"#!{zsh_path}\n")
                fh.writelines(f"#SBATCH --job-name=job_lth{count_}.job\n")
                #fh.writelines(f"#SBATCH --cpus-per-task=30\n")
                fh.writelines(f"#SBATCH --output={job_directory}/job_lth{count_}.out\n")
                # fh.writelines(f"#SBATCH --error=.out/job_lth{count_}.err\n")
                # fh.writelines(f"#SBATCH --workdir={workdir}\n")

                fh.writelines(cmd)

            os.system("sbatch %s" % job_file)
            count_ += 1
        time.sleep(18000)  # 4hours later