import os
from os.path import *
import sys

if __name__ == '__main__':
    args = sys.argv[1:]
    odir = args[0]
    while 1:
        os.system(f'cd {odir}; mcmctree 03_mcmctree.ctl > 03_mcmctree.log')
        if exists(f'{odir}/FigTree.tre'):
            break    
