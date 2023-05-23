"""
The script is used to summary the iqtree outpus and generate in.BV file.
and prepare mcmc_tree.ctl file
"""

import os
from os.path import *
import sys
from glob import glob

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

def get_final_numS(phy,num_s=None):
    if num_s is None:
        num_s = next(open(phy)).split(' ')[0]
    bucket_dict = {}
    c = 0
    for _ in open(phy):
        if _.startswith(num_s):
            c+=1
            bucket_dict[c] = {}
            continue
        else:
            sid = _.split(' ')[0]
            seq = _.split(' ')[-1].strip()
            if len(set(seq))!=1:
                bucket_dict[c][sid] = seq
            # else:
            #     print(set(seq))
    # for k,v in bucket_dict.items():
    #     print(k,len(v))
    return bucket_dict

def process_dir(path):
    if '~' in path:
        return os.path.expanduser(path)
    return path
########################## setting ############################
# software
NW_TOPOLOGY = '/home-user/thliao/anaconda3/bin/nw_topology'
RSCRIPT='/home-user/thliao/anaconda3/bin/Rscript'
RUBY = "/home-user/thliao/bin/ruby"
NW_STATS = '/home-user/thliao/anaconda3/bin/nw_stats'

# template and scripts
REORDER_NODE = "/home-user/thliao/script/sswang_script/reorder_node.rb"
FROM_BS_TO_HESSIAN = '/home-user/sswang/project/Rhizobiales/scripts/dating/hessian/from_bs_to_hessian.R'
template_mcmc = "/home-user/thliao/script/evol_tk/dating_workflow/ctl_template/mcmctree.ctl"
##############################################################

#### params
if __name__ == '__main__':
    args = sys.argv[1:]
    odir,phy_file,calib_tree_file = args
    odir = process_dir(odir)
    phy_file = process_dir(phy_file)
    calib_tree_file = process_dir(calib_tree_file)
    
    ref_tree_file = f"{odir}/ref.tre"
    ####
    def create_inBV(mltree_file, inBV_file, iqtree_outdir):
        no_species = os.popen(f"{NW_STATS} {mltree_file} | grep '^#leaves:' | awk " + " '{print $2}' ").read().strip()
        cmd = f"""echo -e "\n{no_species}\n" > {inBV_file};
        {NW_TOPOLOGY} -bI {mltree_file} >> {inBV_file}
        echo -e "\n" >> {inBV_file}
        cat {iqtree_outdir}/ml.bls >> {inBV_file}
        echo -e "\n" >> {inBV_file}"""
        os.system(cmd)
        gradient = ' '.join(['0'] * (2*int(no_species)-3)) # no. of branches equals 2n-3 where n is the no. of species
        cmd = f"""echo {gradient} >> {inBV_file}
        echo -e "\n" >> {inBV_file}

        echo Hessian >> {inBV_file}
        echo -e "\n" >> {inBV_file}
        cat {iqtree_outdir}/hessian >> {inBV_file} """
        os.system(cmd)

    def single_process(iqtree_outdir):
        """
        convert iqtree output to mcmctree input
        generating in.BV hessian
        """
        boottree_file = join(iqtree_outdir,'iqtree.boottrees')
        mltree_file = join(iqtree_outdir, 'iqtree.treefile')
        cmd = f"{NW_TOPOLOGY} -Ib {boottree_file} | {RUBY} {REORDER_NODE} -i - --ref {ref_tree_file} > {iqtree_outdir}/boot.bls"
        if not exists(f"{iqtree_outdir}/boot.bls"):
            os.system(cmd)
        cmd2 = f"{RUBY} {REORDER_NODE} -i {mltree_file} --ref {ref_tree_file} > {iqtree_outdir}/ml.bls"
        if not exists(f"{iqtree_outdir}/ml.bls"):
            os.system(cmd2)
        cmd3 = f"{RSCRIPT} {FROM_BS_TO_HESSIAN} {iqtree_outdir}/boot.bls {iqtree_outdir}/hessian"
        if not exists(f"{iqtree_outdir}/hessian"):
            os.system(cmd3)
        mltree_file = join(iqtree_outdir, 'iqtree.treefile')
        inBV_file = join(iqtree_outdir,'in.BV')
        create_inBV(mltree_file, inBV_file, iqtree_outdir)

    if not exists(f"{odir}/mcmctree/"):
        os.makedirs(f"{odir}/mcmctree/")
    if not exists(f"{odir}/iqtree/"):
        raise IOError(f'no {odir}/iqtree/ was found')
    for iqtreeodir in glob(odir+'/iqtree/*'):
        single_process(iqtreeodir)
    cmd = f"cat {odir}/iqtree/*/in.BV > {odir}/mcmctree/in.BV"
    os.system(cmd)
    ### copying phy and tree file to mcmctree dir
    cmd = f"cp {phy_file} {odir}/mcmctree/; cp {calib_tree_file} {odir}/mcmctree/"
    os.system(cmd)

    bd_paras = '1 1 0.1'
    rgene_gamma = '1 35 1'
    sigma2_gamma = '1 10 1'
    sampfreq = '30'
    nsample = '20000'
    burnin = str(int(0.1 * (int(sampfreq)*int(nsample))))
    seqfile_b = basename(phy_file)
    treefile_b = basename(calib_tree_file)
    outfile = './03_mcmctree.out'
    seqtype = 2
    clock = 2
    param = {
        'seqfile': seqfile_b,
        'treefile': treefile_b,
        'ndata': str(int(len(get_final_numS(phy_file)))),
        'model':'0',
        'seqtype': seqtype,
        'usedata': "2 in.BV 1",
        'outfile': outfile,
        'clock': clock,
        'BDparas': bd_paras,
        'rgene_gamma': rgene_gamma,
        'sigma2_gamma': sigma2_gamma,
        'burnin': burnin,
        'sampfreq': sampfreq,
        'nsample': nsample,
        'print': '1',
        'alpha': 0.5}

    text = modify(template_mcmc,**param)
    with open(f"{odir}/mcmctree/03_mcmctree.ctl",'w') as f1:
        f1.write(text)
