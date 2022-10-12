#! /env/python

from collections import defaultdict
import os
from tqdm import tqdm
from os.path import *
import click

_t = '/mnt/maple/thliao/data/NCBI/modified_data/annotations/tab.header'
header = [open(_t).read().strip()]

kofamscan_exe = '/home-user/thliao/software/kofamscan/exec_annotation'
ko_list = '/mnt/home-backup/thliao/kofam/20190810/ko_list'

default_kofamscan_exe = '/home-user/thliao/software/kofamscan/exec_annotation'
default_odir= "/mnt/maple/thliao/data/NCBI/modified_data/annotations"
def kofam_anno(faa_list,odir=default_odir,profile=None):
    cmds = []
    for faa in faa_list:
        uid = faa.split('/')[-1].replace('.faa','')
        ofile = f"{odir}/{uid}.kofamout"
        tmp_dir = odir+f'/{uid}_tmp/' 
        if profile is None:
            para = ''
        else:
            para = f" -p {profile} "
        cmd = f"/home-user/thliao/software/kofamscan/exec_annotation {para} --tmp-dir {tmp_dir} -o {odir}/{uid}.kofamout -f mapper-one-line --no-report-unannotated {faa}; rm -rf {tmp_dir}"
        if not exists(ofile):
            cmds.append(cmd)
    return cmds

def correct_reanno(ori_tab,odir,ori_faa=None,tmp_dir='/mnt/ivy/thliao/tmp/',existsed_kos=None,remove=True):
    gid = ori_tab.split('/')[-1].rsplit('.',1)[0]
    if ori_faa is None:
        ori_faa = f"/mnt/maple/thliao/data/NCBI/modified_data/direct_protein_files/{gid}.faa"
    if not exists(ori_tab):
        raise IOError(f'{ori_tab} is not existsed')
    all_kos = [([_ for _ in row.split(' ') if _][2],row)
                for row in open(ori_tab).read().strip().split('\n')
                if not row.startswith('#')]
    ko2rows = defaultdict(list)
    for ko,row in all_kos:
        if existsed_kos is not None:
            if ko not in existsed_kos:
                continue
        ko2rows[ko].append(row)
    tmp_dir = tmp_dir+f'/{gid}'
    if not exists(f"{tmp_dir}/tabular"):
        os.system(f"mkdir -p {tmp_dir}/tabular")
    else:
        os.system(f"rm {tmp_dir}/tabular/*")
    query_faa = realpath(tmp_dir +'/query.faa')
    cmd = f"ln -sf `realpath {ori_faa}` {query_faa} "
    os.system(cmd)
    for k,rows in tqdm(ko2rows.items(),total=len(ko2rows)):
        of = tmp_dir+f'/tabular/{k}'
        with open(of,'w') as f1:
            f1.write('\n'.join(header+rows))
    cmd = f"{kofamscan_exe} -r -k {ko_list} --tmp-dir {tmp_dir} -o {odir}/{gid}.kofamout -f mapper-one-line --no-report-unannotated {query_faa} "
    if remove:
        cmd += f" ; rm -rf {tmp_dir}"
    return cmd


@click.command()
@click.option('-i', 'input_table')
@click.option('-o', 'output_dir')
@click.option('-tmp', 'tmp_dir',default='/mnt/ivy/thliao/tmp/')
@click.option('-faa', 'faa',default=None)
@click.option('-k', 'keep', help='keep the tmp dir or not',default=True, required=False, is_flag=True)
@click.option('-dry_run', 'dry_run', help='run or not',default=False, required=False, is_flag=True)
def main(input_table,output_dir,tmp_dir,keep,dry_run,faa):
    cmd = correct_reanno(input_table,output_dir,
                         ori_faa=faa,tmp_dir=tmp_dir,
                         existsed_kos=None,remove=keep)
    if dry_run:
        with open('./cmds','a') as f1:
            f1.write(cmd+'\n')
    else:
        os.system(cmd)
        
if __name__ == '__main__':
    main()
    # 