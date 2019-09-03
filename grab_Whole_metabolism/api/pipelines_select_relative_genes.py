import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
from Bio import SeqIO
from collections import Counter
import os
from os.path import basename, dirname, abspath,exists,join
from subprocess import check_call


def output_seq_from_df(df, oseq):
    with open(oseq, 'w') as f1:
        for _, row in df.iterrows():
            aa_seq = row['AA seq']
            name = row["locus_name"]
            f1.write(f'>{name}\n')
            f1.write(f"{aa_seq}\n")

def run_cmd(cmd):
    check_call(cmd, shell=True)


def main(infile, target_fa, oseq,project_name):
    df = pd.read_csv(infile, sep='\t')
    output_seq_from_df(df, oseq)
    dbname = abspath(oseq).rpartition('.')[0]
    odir = dirname(oseq)
    if not exists(odir):
        os.makedirs(odir)
    # name
    o1_tab = join(odir,)

    run_cmd(f'diamond makedb --in {oseq} --db {dbname}')
    run_cmd(f"diamond blastp -q {target_fa} -o {o1_tab} -d {dbname} -p 0 -b 5 -c 2")

    tmp_df = pd.read_csv(f'{o1_tab}', sep='\t', header=None)
    records = SeqIO.parse(f'{target_fa}', format='fasta')
    used_gids = set(tmp_df.iloc[:, 0])
    collcect_records = []
    for record in tqdm(records):
        if record.id in used_gids:
            collcect_records.append(record)
    with open(f'{intermedia_faa}', 'w') as f1:
        SeqIO.write(collcect_records, f1, format='fasta-2line')

    run_cmd(f"diamond blastp -q {intermedia_faa} -o {o2_tab} -d /home-user/sswang/db/diamond/kegg/latest/kegg -p 0 -b 5 -c 2")


    pre_df = pd.read_csv(f'{o1_tab}', sep='\t', header=None)
    aft_df = pd.read_csv(f'{o2_tab}', sep='\t', header=None)

    def cal_ratio(locus):
        setA = set(pre_df.loc[pre_df.loc[:, 0] == locus, 1])
        setB = set(aft_df.loc[aft_df.loc[:, 0] == locus, 1])
        intersec = setA.intersection(setB)
        union_set = setA.union(setB)
        if len(intersec) / len(union_set) >= 0.5:
            # real_N_metabolism_genes.append(locus)
            return locus
    locus_set = list(set(pre_df.loc[:, 0]))
    real_N_metabolism_genes = []
    with mp.Pool(processes=64) as tp:
        for locus in tqdm(tp.imap(cal_ratio, locus_set),
                          total=len(locus_set)):
            if locus is not None:
                real_N_metabolism_genes.append(locus)

    real_N_metabolism_genes = set(real_N_metabolism_genes)
    records = SeqIO.parse(f'{intermedia_faa}', format='fasta')
    collect_reads = [_ for _ in records if _.id in set(real_N_metabolism_genes)]

    with open(f'{final_faa}', 'w') as f1:
        SeqIO.write(collect_reads, f1, format='fasta-2line')
