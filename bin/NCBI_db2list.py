import click
from glob import glob
from os.path import *

import kmapper as km
mapper = km.KeplerMapper(verbose=1)

indir = "./"
dirs = glob(f'{indir}/bacteria/*/')
dirs2 = glob(f'{indir}/archaea/*/')


dirs3 = dirs+dirs2
all_ids = [_.split('/')[-2] for _ in dirs3]
with open('./curr_ids','w') as f1:
    f1.write('\n'.join(all_ids))