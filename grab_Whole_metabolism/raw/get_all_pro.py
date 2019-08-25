import os, sys
from glob import glob
from os.path import join, basename, dirname
from tqdm import tqdm
from subprocess import check_call
import pandas as pd
from collections import defaultdict

summary_dict = defaultdict(dict)

indir = "/mnt/home-backup/thliao/metagenomes"

all_fa_gz = glob(join(indir, "**", "*.fa.gz"), recursive=True)
all_fna_gz = glob(join(indir, "**", "*.fna.gz"), recursive=True)

odir = "/mnt/home-backup/thliao/metagenomes/concat_all/"
os.makedirs(odir, exist_ok=True)
for fa in tqdm(all_fa_gz + all_fna_gz):
    ofile = join(odir, 'data', basename(fa).rpartition('.')[0].replace('fna', 'fa'))
    # if not os.path.exists(ofile):
        # check_call(f"gunzip -d -c {fa} > {ofile}", shell=True)
    _list = fa.split('/')
    summary_dict[basename(ofile)]['source'] = _list[_list.index('metagenomes')+1]
all_fa = glob(join(indir, "**", '*.fa'), recursive=True)
all_fna = glob(join(indir, '**', '*.fna'), recursive=True)
total_fa = [_ for _ in (all_fa + all_fna) if 'concat_all' not in _]
for fa in tqdm(total_fa):
    ofile = join(odir, 'data', basename(fa).replace('fna', 'fa'))
    # if not os.path.exists(ofile):
        # check_call(f"cp {fa} {ofile}", shell=True)
    _list = fa.split('/')
    summary_dict[basename(ofile)]['source'] = _list[_list.index('metagenomes')+1]
summary_df = pd.DataFrame.from_dict(summary_dict, orient='index')
summary_df.to_csv(join(odir, 'summary_data.csv'), index=True)

import multiprocessing as mp
def run_cmd(args):
    cmd, odir, sample_name = args
    os.makedirs(f"{odir}/{sample_name}", exist_ok=True)
    if len(os.listdir(f"{odir}/{sample_name}")) != 12:
        check_call(cmd,
                   executable="/home-user/thliao/anaconda3/bin/zsh",
                   shell=True,
                   stdout=open("/dev/null", 'w'),
                   stderr=open("/dev/null", 'w'))
    return 'success %s' % sample_name


prokka_p = "/usr/local/bin/prokka"


def main(idir, odir):
    params = []
    for seq in tqdm(glob(os.path.join(idir, '*.fa'))):
        # sample_name = os.path.basename(seq).split('.')[0]
        sample_name = os.path.basename(seq).rpartition('.')[0]
        cmd = "{prokka_p} {infile} --outdir {odir}/{sample_name} --prefix {sample_name} --force --cpus 5 --quiet"
        params.append((cmd.format(prokka_p=prokka_p,
                                  infile=seq,
                                  odir=odir,
                                  sample_name=sample_name, ),
                       odir,
                       sample_name))
    with mp.Pool(processes=13) as tp:
        for r in tqdm(tp.imap(run_cmd, params),
                      total=len(params)):
            t = r


main("/mnt/home-backup/thliao/metagenomes/concat_all/data",
     "/mnt/home-backup/thliao/metagenomes/concat_all/prokka_o")
