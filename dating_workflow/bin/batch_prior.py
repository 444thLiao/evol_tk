import os
from glob import glob
from os.path import *

import click

from dating_workflow.bin.dating_pro import run_nodata_prior, process_path, get_num_phy_file, run, tqdm, mp


@click.command()
@click.option('-i', '--in_phy', 'in_phyfile',help="single alignment file with phylo format")
@click.option('-it', '--in_tree', 'in_treefile',help="single tree file with newick format with calibration information")
@click.option('-o', 'odir',help="output directory")
@click.option('-nucl', 'use_nucl', is_flag=True, default=False,help="If you use nucleotide sequence, you should pass is. But if you pass it when you use amino acids, if will raise silent errors.")
@click.option('-sf', 'sampfreq', default='2')
@click.option('-p', 'print_f', default='2')
@click.option('-rg', 'rgene_gamma', default='1 35 1')
@click.option('-sg', 'sigma2_gamma', default='1 10 1')
@click.option('-c', 'clock', default='2')
@click.option('-np', 'num_parellel', default=10, help="Number of processes you want to parallel")
@click.option('-f', 'force', is_flag=True, default=False,
              help="Overwrite previous results or not.")
def cli(in_phyfile, in_treefile,
        odir,
        use_nucl,
        num_parellel,
        sampfreq, print_f, rgene_gamma, sigma2_gamma, clock,force):
    in_phyfile = process_path(in_phyfile)
    ndata = get_num_phy_file(in_phyfile)
    params_dict = {'sampfreq': str(sampfreq),
                   'print': str(print_f),
                   'rgene_gamma': rgene_gamma,
                   'sigma2_gamma': sigma2_gamma,
                   'clock': clock}

    if '*' in in_treefile:
        trees = [process_path(_) for _ in glob(in_treefile)]
    else:
        trees = [process_path(in_treefile)]
    cmds = []
    for tree in trees:
        name = basename(tree).rpartition('.')[0]
        new_odir = join(odir, name)
        if not exists(new_odir):
            os.makedirs(new_odir)
        cmd = run_nodata_prior(in_phyfile=in_phyfile,
                               in_treefile=tree,
                               odir=new_odir,
                               ndata=ndata,
                               use_nucl=use_nucl,
                               params_dict=params_dict)
        if (not exists(join(new_odir,'FigTree.tre'))) or force:
            cmds.append(cmd)
    with mp.Pool(processes=num_parellel) as tp:
        _ = list(tqdm((tp.imap(run, cmds)), total=len(cmds)))


if __name__ == '__main__':
    cli()

    # python3 ~/script/evolution_relative/dating_workflow/bin/batch_prior.py -i ./dating_for/phy_files/83g_nuc_concat.phy -it './dating_for/cal_tree/83g_set*.newick' -o ./dating_for/83g/batch_prior_nucl -nucl -sf 20 -p 1 -rg '1 30 1' -c 2
