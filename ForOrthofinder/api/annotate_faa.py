import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from subprocess import check_call
import multiprocessing as mp
import click
from glob import glob
import os
from tqdm import tqdm
import pandas as pd
import io
from hmmparser import HMMparser
from collections import Counter
from pandas.errors import EmptyDataError

default_db = "/home-db/pub/protein_db/TIGRFAM/14.0_release/TIGRFAM.HMM"


def run_cmd(seqfile, db, tbloutput, evalue, cpu):
    if os.path.exists(tbloutput):
        pass
    else:
        template = f"hmmscan --tblout {tbloutput} --acc --noali --notextw -E {evalue} --cpu {cpu} {db} {seqfile} "
        # print(template)
        check_call(template,
                   shell=True,
                   stdout=open("/dev/null", 'w'))
    try:
        hmm = HMMparser(tbloutput)
        txt_table = ''
        for row in hmm.matrix:
            if len(row) > 22:
                row = row[:22]
            elif len(row) < 22:
                row = row + [''] * (22 - len(row))
            txt_table += '\t'.join(row) + '\n'
        with open(tbloutput, 'w') as f1:
            print(txt_table, file=f1)
    except KeyError:
        # it means that it have been processed yet.
        txt_table = open(tbloutput, 'r').read()
    try:
        table = pd.read_csv(io.StringIO(txt_table),
                            sep='\t',
                            header=None,
                            index_col=None)
    except EmptyDataError:
        return pd.DataFrame(), seqfile
    return table, seqfile


def exec_cmd(args):
    return run_cmd(*args)


def get_seqs(indir, suffix):
    seqs = glob(os.path.join(indir, "*." + suffix))
    if not seqs:
        print(os.path.join(indir, "*." + suffix), "doesn't match any files")
    return seqs


def process_table(table: pd.DataFrame):
    if table.shape[0] == 0:
        return '?', '?'
    sorteed_table = table.sort_values(4)  # evalue for full sequence
    droppedDup_table = sorteed_table.drop_duplicates(2)
    collect_id = droppedDup_table.loc[:, 1]
    count_ = Counter(collect_id)
    most_id = list(sorted(count_, key=lambda x: count_[x]))[-1]
    _cache = table.loc[table.loc[:, 1] == most_id, 18]
    if _cache.shape[0] !=0:
        gene_name = list(_cache)[0].strip(':')
    else:
        gene_name = "?"
    return most_id, gene_name


@click.command()
@click.option("-i", "inpath", help="single file or a directory")
@click.option("-s", "suffix", help="suffix when -i is a directory", default='faa')
@click.option("-o", "output", default=None, help="optional")
@click.option("-n", "cpu", default=2)
@click.option("-p", "parallel", default=4)
@click.option("-db", default=default_db)
@click.option("-e", "evalue", default="1e-10")
def cli(inpath, suffix, output, cpu, db, evalue, parallel):
    inpath = os.path.abspath(inpath)
    if os.path.isdir(inpath):
        seqs = get_seqs(inpath, suffix)
    else:
        seqs = [inpath]
    args = []
    for seq in seqs:
        if output is None:
            output_p = seq.replace('.%s' % suffix,
                                   ".hmmscan")
        elif os.path.isdir(output):
            output_p_base = os.path.basename(seq).replace('.%s' % suffix,
                                                          ".hmmscan")
            output_p = os.path.join(output, output_p_base)
        else:
            output_p = output
        args.append((seq,
                     db,
                     output_p,
                     evalue,
                     cpu))
    f1 = open(os.path.join(os.path.dirname(output_p), 'summary_hmm.log'), 'w')
    print("\t".join(["OG_id",
                     "num_genomes",
                     "annotated_id",
                     "annotated_gene_name"]),
          file=f1)
    with mp.Pool(processes=parallel) as tp:
        for table, seqfile in tqdm(tp.imap(exec_cmd, args)):
            most_id, gene_name = process_table(table)
            num_genomes = table.drop_duplicates(2).shape[0]
            print('\t'.join([os.path.basename(seqfile),
                             str(num_genomes),
                             most_id,
                             gene_name]),
                  file=f1)
            f1.flush()
    f1.close()


if __name__ == '__main__':
    cli()
