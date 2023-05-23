

import os
from os.path import *
import sys
sys.path.insert(0, '/home-user/thliao/script/evol_tk')
from glob import glob
from bin.multiple_sbatch import sbatch_all

def out_phy(outfile,p2s2seq):
    with open(outfile,'w') as f1:
        for v,s2seq in p2s2seq.items():
            num_s = len(s2seq)
            num_sites = len(list(s2seq.values())[0])
            f1.write(f"{num_s}     {num_sites}\n")
            for k,v in s2seq.items():
                f1.write(f"{k}            {v}\n")

def get_final_numS(phy,num_s=None):
    if num_s is None:
        num_s = next(open(phy)).strip(' \n').split(' ')[0]
    bucket_dict = {}
    c = 0
    for row in open(phy):
        row = row.strip('\n')
        if not row:continue
        if row.startswith(num_s):
            c+=1
            bucket_dict[c] = {}
            continue
        else:
            sid = row.split(' ')[0]
            seq = row.split(' ')[-1].strip()
            if len(set(seq))!=1:
                bucket_dict[c][sid] = seq
    return bucket_dict

def process_dir(path):
    if '~' in path:
        return os.path.expanduser(path)
    return path
################# default setting #################
IQTREE = '/home-user/thliao/software/iqtree-2.1.3-Linux/bin/iqtree2'
cpu = 8
_extra_argu = ''  # will be passed to iqtree2
_bootstrap_argu = '-b 1000'
###################################################


def main(model,phy_file,odir,in_BV,dry_run=False):
    # model = 'LG+G+C60'
    # model_name = 'C60'
    # phy_file = ""
    # odir = ''
    # in_BV = ''  
    cmds = []
    # testing and creating output dir
    if not exists(phy_file):
        raise IOError(f"no {phy_file}")
    if not exists(odir):
        os.makedirs(odir)
    if not exists(in_BV):
        raise IOError(f"no {in_BV}")
    model_argu = f'-m {model}'
    ## step1: generating iqtree commands
    cmd = f"sed '4!d' {in_BV} > {odir}/ref.tre"
    os.system(cmd)
    ref_tree_file = f"{odir}/ref.tre" 
    te_argu = f'-te {ref_tree_file}'
    b_d = get_final_numS(phy_file)
    iqtree_outdir = odir+f'/iqtree/'
    if not exists(iqtree_outdir):
        os.makedirs(iqtree_outdir)
        os.makedirs(odir+f'/mcmctree')
        
    for k,s2seq in b_d.items():
        ali_file1 = f"{iqtree_outdir}/{k}/aln.phy"
        if not exists(dirname(ali_file1)):
            os.makedirs(dirname(ali_file1))
        out_phy(ali_file1,{'0':s2seq}) # output phy file
        if not model.split('+')[-1].startswith('C'):
            # not using C60 model
            cmd = f"""{IQTREE} -redo -s {ali_file1} -pre {iqtree_outdir}/{k}/iqtree -nt {cpu} -quiet {model_argu} {_bootstrap_argu} {te_argu} {_extra_argu}"""
        else:
            # use C60 model
            cmd = f"""{IQTREE} -redo -s {ali_file1} -pre {iqtree_outdir}/{k}/guide -nt {cpu} -quiet -m LG4M+G {te_argu} {_extra_argu} ; {IQTREE} -redo -s {ali_file1} -pre {iqtree_outdir}/{k}/iqtree -nt {cpu} -quiet {model_argu} {_bootstrap_argu} {te_argu} {_extra_argu} -ft {iqtree_outdir}/{k}/guide.treefile"""
        cmds.append(cmd)
    if not dry_run:
        sbatch_all(cmds,thread_per_tasks=cpu,job_name='IQtree')
    else:
        with open(f"{odir}/iqtree.cmd",'w') as f1:
            f1.write('\n'.join(cmds))
    

    # calib_tree_file = strategies_params[setname][1]
    # cmd = f"python /home-user/thliao/script/dating_raw/hessian_s2_bin.py {odir} {phy_file} {calib_tree_file}"
    # cmds.append(cmd)

    # cmd = f"python /home-user/thliao/script/evol_tk/dating_workflow/bin/dating_repeat.py -i {odir+f'/mcmctree'} -o {odir+f'/CUSTROM_model_MCMC'} -it {calib_tree_file}"
    # cmds.append(cmd)
    # dating_repeat.py will use the previous the same configuration as the previous mcmc_output 
   
# in_BV = ''
# phy_file = "/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/refined_genes_DeltaL/sliding_windows/P39_10/move10_237.57/RP39_pf/3pf.phy"
# model = 'LG+G+C60'
# odir = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/refined_genes_DeltaL/sliding_windows/P39_10/move10_237.57/testC60'
# calib_tree_file = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/system_dating/cal/P39_B5P3.newick'
    
import click
@click.command()
@click.option('--model', default='LG+G+C60', help='model name')
@click.option('--phy_file', help='phylip formatted file. Can be multiple partitions')
@click.option('--odir', help='output dir. It will automatically generate a subfolder named "mcmctree" and subfolder named "iqtree"')
@click.option('--inbv', help='a pre-calculated in.BV file. It mainly used for extracting a tree file. Thus, it should share the same topology as your used tree file finally. ')
@click.option('--dryrun',is_flag=True, default=False, help='directly run it or generate a file containing all commands.')
def cli(model,phy_file,odir,inbv,dryrun):
    odir = process_dir(odir)
    phy_file = process_dir(phy_file)
    inbv = process_dir(inbv)
    main(model,phy_file,odir,inbv,dryrun)


if __name__ == '__main__':
    cli()
    # usages:
    # python customModel_MCMC.py --model 'LG+G+C60' --phy_file '/home-user/thliao/script/dating_raw/example_data/3pf.phy' --odir '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/refined_genes_DeltaL/sliding_windows/P39_10/move10_237.57/testC60' --inbv '/home-user/thliao/script/dating_raw/example_data/in.BV' --dryrun    
