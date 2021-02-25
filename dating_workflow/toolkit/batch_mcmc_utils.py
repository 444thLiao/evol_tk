import os
import pandas as pd
from os.path import abspath, basename, dirname, join,exists
from api_tools.for_tree.format_tree import add_cal_api
from glob import glob
from dating_workflow.bin.dating_pro import modify, run

def generate_cal(in_xlsx,
                 odir,
                 colname2groups,
                 colname2name={}
                 ):
    calinfodf = pd.read_excel(in_xlsx,index_col=0)
    cal_files = []
    cal_names = []
    for cal_name, row in calinfodf.iterrows():
        rows = []
        for node, t in row.items():
            if not pd.isna(t):
                rows.append(f"{colname2groups[node]}\t{str(t)}\t{colname2name[node]}")
        cal_file =join(odir,f'cal_{cal_name}.txt') 
        with open(cal_file, 'w') as f1:
            f1.write('\n'.join(rows))
        cal_files.append(cal_file)
        cal_names.append(cal_name)
    return list(zip(cal_names,cal_files))


def generate_cal_trees(cal_sets,
                       in_tree,
                       odir,
                       format=3,
                       redo = False,
                       extra_suffix=''):

    """generate a set of cal_tree with an single input tree and batch of calibration files.
    
    Args:
        cal_sets ([list]): list of zipped (cal_names,cal_files)
        in_tree ([str]): [description]
        format (int, optional): format for reading the in_tree, see format descriptions of ete3. Defaults to 3.
        redo (bool, optional): redo or not for existing calibration tree. If redo, it will overlap the exised tree file. Defaults to False.
        extra_suffix (str, optional): it would be added to at the front of the based name of calibrated tree and separated with '_' to the name of calibration set  . Defaults to ''.

    Returns:
        [dict]: dict of trees with calibration information
    """
    cal_trees = {}
    for (cal_name,cal_file) in cal_sets:
        otree = join(odir,f"{extra_suffix+'_'}{cal_name}")
        if redo or not exists(otree):
            add_cal_api(in_tree,
                        otree,
                        cal_file,
                        format=format)
            cal_trees[cal_name] = otree
    return cal_trees


def generate_batch_mcmc(cal_trees,
                        template_dir,
                        odir,
                        program_name='',
                        clock_t= '2',
                        seqtype='nucl',
                        seqfile=None,
                        print=1,
                        sampfreq=20,
                        nsample=20000,
                        burnin=2000
                        ):
    """generate batch mcmc directory according to the template dir which contains the precalculated in.BV file.

    Args:
        cal_trees (dict): dict of calibrated trees with the set name as their key. 
        template_dir (str): dir with in.BV
        odir (str): ouput directroy, it will automated generate descending directory following f"{odir}/{program_name}/{seqtype}/clock{clock_t}". 
        program_name (str, optional): for generated ouput directory. Defaults to ''.
        clock_t (str, optional): 3 for AR. 2 for IR . Defaults to '2'.
        seqtype (str, optional): [description]. Defaults to 'nucl'.
        seqfile (str, optional): [description]. Defaults to None.
        print (int, optional): [description]. Defaults to 1.
        sampfreq (int, optional): [description]. Defaults to 20.
        nsample (int, optional): [description]. Defaults to 20000.
        burnin (int, optional): [description]. Defaults to 2000.

    Returns:
        list: list of commands suitable to run mcmc
    """    
    odir = f"{odir}/{program_name}/{seqtype}/clock{clock_t}"
    cmds = []
    for cal_name,cal_tree in cal_trees.items():
        set_name = cal_name
        pre_ctl = glob(f'{template_dir}/*.ctl')
        if not pre_ctl:
            continue
        else:
            pre_ctl = pre_ctl[0]
        inbv = abspath(glob(f'{template_dir}/in.BV')[0])
        param = {
            # 'seqfile': seqfile_b,
            'treefile': abspath(cal_tree),
            # 'ndata': ndata,
            # 'seqtype': seqtype,
            'usedata': "2 in.BV 1",
            'outfile': './03_mcmctree.out',
            'clock': clock_t,
            # 'BDparas': bd_paras,
            'rgene_gamma': '1 100 1' if seqtype=='nucl' else '1 30 1',
            # 'sigma2_gamma': sigma2_gamma,
            'burnin': burnin,
            'sampfreq': sampfreq,
            'nsample': nsample,
            # 'alpha': 0.5,
            'print': print
        }
        if seqfile is not None:
            param['seqfile'] = seqfile
            
        # modify these ctl
        for repeat_n in ['run1','run2']:
            onew_name = f'{set_name}_{repeat_n}'
            if exists(pre_ctl):
                # prepare dir and soft link the in.BV
                final_odir = join(odir, onew_name)
                os.makedirs(final_odir, exist_ok=True)
                os.system(f"ln -sf {inbv} {final_odir}/ ")
                output_ctl = join(final_odir, '04_mcmctree.ctl')
                # generate the text
                text = modify(pre_ctl, **param)
                with open(output_ctl, 'w') as f1:
                    f1.write(text)
                
                # generate the command
                cmd = f"cd {join(final_odir)} ; mcmctree ./04_mcmctree.ctl > ./run.log "

                # check the completeness of previous run
                if not exists(join(final_odir, 'FigTree.tre')):
                    cmds.append(cmd)
                    if exists(join(final_odir, 'mcmc.txt')):
                        # if not figTree.tre but mcmc.txt, it means it doesn't complete at last time
                        os.system(f"rm {join(final_odir, 'mcmc.txt')} ")
    return cmds
