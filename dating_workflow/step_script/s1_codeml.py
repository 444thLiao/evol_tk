
from ete3 import Tree
import click
from glob import glob
import multiprocessing as mp
from subprocess import check_call,check_output
import subprocess
import os
from os.path import *
from tqdm import tqdm
indir = '/home-user/thliao/data/nitrification_for/dating_for/mcmc_for/proteins'
in_phyfile = './concat_aln.phy'
in_treefile = './243g_120gene.calibrations.newick'


# in_phyfile = './concat_aln.phy'
# in_treefile = './iqtree_sorted_topology.newick'

# def separate_phy(in_phyfile):
#     part_collect = []
#     _cache = []
#     for row in open(in_phyfile):
#         if not row.strip('\n'):
#             continue
#         if not row.startswith('GCA') and row.split(' ')[0].isnumeric():
#             if _cache:
#                 part_collect.append(_cache)
#             _cache = []
#         _cache.append(row.strip('\n'))
#     if _cache:
#         part_collect.append(_cache)
#     parts = ['\n'.join(_) for _ in part_collect]
#     return parts

def modify(file,**kwargs):
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
    if isinstance(args,str):
        cmd = args
    else:
        cmd,log = args
    try:
        subprocess.check_output(cmd,shell=1)
    except subprocess.CalledProcessError as e:
        print('error')
        error = e.output
        return error
        
    # if isinstance(args,str):
    #     cmd = args
    #     check_call(cmd,shell=1,)
    #            #stdout=open('/dev/null','w'))
    # else:
    #     cmd,log = args
    #     check_call(cmd,shell=1,
    #            stdout=open(log,'w'),
    #            stderr=open(log,'w'))
    

# parts = separate_phy(in_phyfile)
template_ctl_01 = './01_mcmctree.ctl'
new_01_ctl = './01_mcmctree_modify.ctl'
params = {'ndata':25,
          'seqfile':in_phyfile,
          'treefile':in_treefile,
          'outfile':'01_out'}
text = modify(template_ctl_01,**params)
with open(new_01_ctl,'w') as f1:
    f1.write(text)
run(f"export PATH=''; /home-user/thliao/software/paml4.9j/bin/mcmctree {new_01_ctl}")

ctls = glob('./tmp0*.ctl')
for ctl in ctls:
    name = basename(ctl).replace('.ctl','')
    os.makedirs(f'./{name}',exist_ok=1)
    os.system(f'mv {name}.* {name}/')
    

params = []
ctls = glob('./*/tmp0*.ctl')
ctls = [_ for _ in ctls if 'modify' not in _]
for ctl in ctls:
    new_text = modify(ctl,
                      **{'model':2,'aaRatefile':'../lg.dat','fix_alpha':0,'alpha':0.5,'ncatG':4})
    new_file = ctl.replace('.ctl','.modify.ctl')
    with open(new_file,'w') as f1:
        f1.write(new_text)
        params.append(f"cd {dirname(new_file)}; /home-user/thliao/software/paml4.9j/bin/codeml {basename(new_file)} > {basename(new_file).replace('.modify.ctl','.log')}")

with mp.Pool(processes= 30) as tp:
    r = list(tqdm((tp.imap(run,params)),total=len(params)))


text = ''
for rst2 in sorted(glob("./tmp*/rst2"),
                   key=lambda x: int(x.split('/')[1].replace('tmp',''))):
    text += open(rst2).read()
with open('./in.BV','w') as f1:
    f1.write(text)
# run("cat ./tmp*/rst2 > in.BV")


# param = {}
# text = modify('./01_mcmctree.ctl',
#                **param)




# for final mcmctree
bd_paras = '1 1 0.1'
rgene_gamma = '1 35 1'
sigma2_gamma = '1 10 1'
burnin = '2000'
sampfreq = '2'
nsample = '20000'
seqfile_b = '.'+in_phyfile
treefile_b = '.'+in_treefile
outfile = './out'
ndata = 23
seqtype = 2
clock = 2
param = {'seqfile':seqfile_b, 
         'treefile':treefile_b, 
         'ndata':ndata, 
         'seqtype':seqtype, 
         'usedata':"2 in.BV 1", 
         'outfile':outfile,
         'clock':clock, 
         'BDparas':bd_paras, 
         'rgene_gamma':rgene_gamma, 
         'sigma2_gamma':sigma2_gamma, 
         'burnin':burnin, 
         'sampfreq':sampfreq, 
         'nsample':nsample, 
         'alpha':0.5}
text = modify('./01_mcmctree.ctl',
               **param)
os.makedirs('./03step',exist_ok=1)
os.system('cp ./in.BV ./03step/')
ofile = './03step/03_mcmctree.ctl'
with open(ofile,'w') as f1:
    f1.write(text)
run(f"/home-user/thliao/software/paml4.9j/bin/mcmctree {ofile} > {ofile.replace('.ctl','.log')}")


def cli(tree,aln_dir,odir):
    pass

if __name__ == "__main__":
    cli()