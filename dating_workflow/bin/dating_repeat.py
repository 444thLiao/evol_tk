
import os
from glob import glob
from os.path import *

import click
from tqdm import tqdm
from dating_workflow.bin.dating_pro import modify


def repeat_mcmc(indir,intree,odir,params_dict={},same=True):
    if not exists(odir):
        os.makedirs(odir)
    ctl = glob(join(indir,'*.ctl'))[0]
    phy = glob(join(indir,'*.phy'))[0]
    cmd = f"cp {ctl} {phy} in.BV {intree} {odir}/"
    os.system(cmd)
    ofile = join(odir,basename(ctl))
    
    if same:
        params_dict = {}
    if params_dict:
        text = modify(ctl,
                      **params_dict)
        ofile = join(odir, '03_mcmctree.ctl')
        with open(ofile, 'w') as f1:
            f1.write(text)
    tqdm.write("start running the final mcmctree. ")            
    params = [(f"cd {odir}; mcmctree 03_mcmctree.ctl 2>&1",
              ofile.replace('.ctl', '.log'))]
    while 1:
        os.system(params[0][0]+' > ' +ofile.replace('.ctl', '.log'))
        if exists(join(dirname(ofile),'FigTree.tre')):
            break


@click.command()
@click.option('-i', '--indir', 'indir',help='mcmc_for directory')
@click.option('-it', '--intree', 'in_treefile',help='newick format tree with calibration information')
@click.option('-o', 'odir',help='')
@click.option('-s', 'same', is_flag=True, default=True,help='same as previous')
@click.option('-sf', 'sampfreq', default='20',help="sample frequency  [20]")
@click.option('-p', 'print_f', default='2',help="verbose of print  [2]")
@click.option('-rg', 'rgene_gamma', default='1 35 1',help="rgene_gamma: prior on mutation rate   [1 35 1]")
@click.option('-sg', 'sigma2_gamma', default='1 10 1',help="sigma2_gamma: shape and scale parameters  [1 10 1]")
@click.option('-bd', 'bdparse', default='1 1 0.1',help="verbose of print  [2]")
@click.option('-c', 'clock', default='2',help="2 indicate using IR clock model, while 3 denote AR clock model")
def cli(indir, in_treefile, 
        odir, same,sampfreq,
        print_f, rgene_gamma, sigma2_gamma,
        bdparse,
        clock):
    params_dict = {'sampfreq': str(sampfreq),
                   'print': str(print_f),
                   'rgene_gamma': rgene_gamma,
                   'sigma2_gamma': sigma2_gamma,
                   "BDparas":bdparse,
                   'clock': clock}
    repeat_mcmc(indir,in_treefile,odir,params_dict,same)


if __name__ == "__main__":
    cli()



