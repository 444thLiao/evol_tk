import sys
import os

import time
def get_tasks(f):
    lines = open(f).read().split('\n')
    lines = [_ for _ in lines if _]
    os.system('mkdir -p ./tmp/')

    used_files = []
    for cmd in lines:
        curr_time = str(int(time.time()))
        sh_file = f'./tmp/{curr_time}.sh'
        with open(sh_file,'w') as f1:
            f1.write("#!/bin/zsh\n")
            f1.write(cmd)
        os.system(f'sbatch {sh_file}')
        used_files.append(sh_file)
        # os.system(f'rm {sh_file}')

if __name__ == '__main__':
    argv = sys.argv[1]