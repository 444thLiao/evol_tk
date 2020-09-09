#!/usr/bin/env python

import os
import sys
def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        print(dir)
        os.mkdir(dir)

if __name__ == "__main__":
    job_directory = f"{os.getcwd()}/.job" 
    # Make top level directories
    mkdir_p(job_directory)

    cmd_files = sys.argv[1]
    cmds=[_ for _ in open(cmd_files).read().split('\n')]

    count_ = 0
    for cmd in cmds:
        workdir = cmd.split(';')[0].strip().split(' ')[-1]
        cmd = cmd.split(';')[-1]
        job_file = os.path.join(job_directory,f"job_lth{count_}.job" )


        with open(job_file,'w') as fh:
            fh.writelines("#!/home-user/thliao/anaconda3/bin/zsh\n")
            fh.writelines(f"#SBATCH --job-name=job_lth{count_}.job\n")
            # fh.writelines(f"#SBATCH --output=.out/job_lth{count_}.out\n")
            # fh.writelines(f"#SBATCH --error=.out/job_lth{count_}.err\n")
            fh.writelines(f"#SBATCH --workdir={workdir}\n")
            
            fh.writelines(cmd)

        os.system("sbatch %s" %job_file)
        count_ += 1
    