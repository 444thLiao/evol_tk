import multiprocessing as mp
import os
import subprocess
from glob import glob
from os.path import *
from subprocess import check_call

import click
from tqdm import tqdm

template_dir = join(dirname(dirname(__file__)), 'ctl_template')
codeml_ctl = join(template_dir, 'codeml.ctl')
aaRatefile = join(template_dir, 'lg.dat')

paml_bin = "/home-user/thliao/software/paml4.9j/bin"




def modify(file, **kwargs):
    text = open(file).read()
    text = text.split('\n')
    new_text = []
    for row in text:
        key = row.split('=')[0].strip()
        if key in kwargs:
            new_text.append(f"{key} = {kwargs[key]}")
        else:
            new_text.append(row)
    return '\n'.join(new_text)


def run(args):
    if isinstance(args, str):
        cmd = args
        log = '/dev/null'
    else:
        cmd, log = args
    try:
        # subprocess.check_output(cmd, shell=1)
        check_call(cmd,
                   shell=1,
                   stdout=open(log, 'w'))

    except subprocess.CalledProcessError as e:
        print('error', e.output)
    if log != '/dev/null':
        t = open(log, 'r', newline='\n').read().replace('\r', '\n')
        with open(log, 'w') as f1:
            f1.write(t)
            
odir = './SubRate'
phy_files = '195g_concat.phy'
gids = ''
