
from ete3 import Tree
import click
from glob import glob
import multiprocessing as mp
from subprocess import check_call
import os
from os.path import *

in_phyfile = '/home-user/thliao/data/nitrification_for/dating_for/bac120_annoate/concat/255g/concat_aln.phy'
in_treefile = '/home-user/thliao/data/nitrification_for/dating_for/bac120_annoate/concat/255g/iqtree.treefile'




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

# param = {}
# text = modify('./01_mcmctree.ctl',
#                **param)

def run(cmd):
    check_call(cmd,shell=1)
params = []
ctls = glob('./tmp*.ctl')
for ctl in ctls:
    if 'modify.' in ctl:
        continue
    new_text = modify(ctl,
                      **{'model':2,'aaRatefile':'../lg.dat','fix_alpha':0,'alpha':0.5,'ncatG':4})
    new_file = join(ctl.replace('.ctl',''),
                    basename(ctl))
    os.makedirs(ctl.replace('.ctl',''),exist_ok=1)
    run(f"cp {ctl.replace('ctl','*')} {ctl.replace('.ctl','')}/")
    with open(new_file,'w') as f1:
        f1.write(new_text)
        params.append('cd %s; /home-user/thliao/software/paml4.9j/bin/codeml %s ' % (ctl.replace('.ctl',''),
                                                                                     basename(ctl)))

    
with mp.Pool(processes= 30) as tp:
    list(tp.imap(run,params))

run("cat tmp*/rst2 > in.BV")


# for final mcmctree
bd_paras = '1 1 0.1'
rgene_gamma = '1 35 1'
sigma2_gamma = '1 10 1'
burnin = '2000'
sampfreq = '2'
nsample = '20000'
seqfile_b = '../concat_aln.phy'
treefile_b = '../iqtree_sorted_topology.newick'
ndata = 23
seqtype = 2
clock = 2
param = {'seqfile':seqfile_b, 
         'treefile':treefile_b, 
         'ndata':ndata, 
         'seqtype':seqtype, 
         'usedata':"2 in.BV 1", 
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

with open('./03_mcmctree.ctl','w') as f1:
    f1.write(text)
run('/home-user/thliao/software/paml4.9j/bin/mcmctree %s' % './03_mcmctree.ctl')


def cli(tree,aln_dir,odir):
    pass

if __name__ == "__main__":
    cli()