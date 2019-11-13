#################################################################################
####  Extract protein sequences from output of orthofinder
####  This script is mainly focus on the low number of genomes found at single genoms
####  It mainly because input genomes shared little bit of genes (single cell sequencing? virus evolution?)
####
#################################################################################
import sys
import os
from os.path import dirname

sys.path.insert(0, dirname(dirname(__file__)))
import pandas as pd
from Bio import SeqIO
from os.path import *
import click
from toolkit.utils import get_dict, get_protein, get_single_copy, get_summary_statistic
from tqdm import tqdm
from collections import Counter
from api.getSeqofOG import get_seq_with_OG
from glob import glob
from toolkit.utils import run_cmd


def choose_paralog():
    pass


def detect_too_lowgroup(num_each_g):
    total = sum(num_each_g)
    too_low = [(k, v)
               for k, v in num_each_g.items()
               if v <= total * 0.01]
    if too_low:
        print("detect too low group %s, it will be removed" % str(too_low))
    return [_[0] for _ in too_low]


def select_OG(data, rr, rn,rgf, total, group_info):

    if total <= 0:
        raise Exception("Wrong input total number of genomes")
    # selected genomes list
    if rgf:
        selected_genomes = open(rgf).read().split('\n')
        selected_genomes = [_
                            for _ in selected_genomes
                            if _ and _ in data.columns]
        
        data = data.loc[:,selected_genomes]
        total = len(selected_genomes)
    # remained ratio/numbers, finally convert to number
    if rr is None and rn is None:
        
        raise Exception("no input params")
    if rr:
        rn = int(total * rr)
    else:
        rr = float(rr / total)
    print("with input parameters, extract OG which present at %s genomes" % rn)

    if group_info is None:
        # if group info ,it will count by number separated
        count_series = data.count(1)
        keep_OG = list(count_series.index[count_series >= rn])
        selected_genomes = data.loc[keep_OG, :]
        selected_genomes = list(selected_genomes.loc[:, ~selected_genomes.isna().all(0)].columns)
        return keep_OG, selected_genomes
    else:
        num_each_g = group_info.value_counts()
        num_each_g = num_each_g.drop(detect_too_lowgroup(num_each_g))
        # drop too low manually assigned group
        id2g = group_info.to_dict()
        remained_num_each_g = (num_each_g * rr).astype(int).to_dict()

        keep_OG = []
        copy_data = data.copy()
        copy_data.columns = list(map(lambda x: id2g.get(x, None),
                                        list(copy_data.columns)))
        for rid, row in copy_data.iterrows():
            subset_row = row[~row.isna()]
            count_subset = Counter(subset_row.index)
            # print(count_subset)
            if all([count_subset.get(g, 0) >= v
                    for g, v in remained_num_each_g.items()]):
                keep_OG.append(rid)
        print("Found required %s OG " % len(keep_OG))

        selected_subset = data.loc[keep_OG, :]
        selected_subset = selected_subset.loc[:, ~selected_subset.isna().all(0)]
        return keep_OG, list(selected_subset.columns)

def get_group_info(group_file,
                   group_column, ):
    if group_file is None and group_column is None:
        return
    data = pd.read_csv(group_file, sep=',', index_col=0, low_memory=False)
    if '\t' in data.columns[0]:
        data = pd.read_csv(group_file, sep=',', index_col=0, low_memory=False)
    if group_column not in data.columns:
        raise Exception("Wrong column provided '%s'" % group_column)
    group_info = data.loc[:, group_column]
    return group_info


def do_mafft(indir, suffix='faa'):
    for p_file in tqdm(glob(join(indir, '*.' + suffix))):
        pre_name = p_file.replace('.%s' % suffix,
                                  '')
        if (not exists(f"{pre_name}.aln")) or (getsize(f"{pre_name}.aln") ==0):
            run_cmd(f"mafft --maxiterate 1000 --genafpair --thread -1 {pre_name}.{suffix} > {pre_name}.aln")


@click.command()
@click.option("-i", "infile",help="normally is the directory called `Results_Nov09_1/Orthogroups/Orthogroups.tsv`")
@click.option("-o", "output_dir", help="the directory of output to")
@click.option("--only_ortholog", help="use ortholog only or use paralog? default only use ortholog", is_flag=True, default=True)
@click.option("-rr", "remained_ratio", help='choose ', default=1.0,type=float)
@click.option("-rn", "remained_num", help='choose ', default=None)
@click.option("-rgf", "remained_genomes_file", help="selecting single copy OG according a file which each line is a genome", default=None)
@click.option("-g", "group_file", help='choose ', default=None)
@click.option("-gc", "group_column", default=None)
@click.option("-doMSA", "doMSA", help='do MSA with mafft or not?', default=False, is_flag=True)
def main(infile,
         output_dir,
         only_ortholog,
         remained_ratio,
         remained_num,
         remained_genomes_file,
         group_file,
         group_column,
         doMSA
         ):
    if doMSA:
        tqdm.write("performing the MSA process")
        do_mafft(output_dir)
        return 
    if only_ortholog:
        data = get_single_copy(infile)
        group_info = get_group_info(group_file, group_column)
        get_OGs, result_genome = select_OG(data,
                                           rr=remained_ratio,
                                           rn=remained_num,
                                           rgf=remained_genomes_file,
                                           total=data.shape[1],
                                           group_info=group_info,
                                           )
        os.makedirs(output_dir, exist_ok=True)
        with open(join(output_dir, 'selected_genomes.txt'), 'w') as f1:
            f1.write('\n'.join(result_genome))
        tqdm.write("Extracting required %s OG" % len(get_OGs))
        get_seq_with_OG(infile,
                        get_OGs,
                        genomes_list=result_genome,
                        output_dir=output_dir,
                        single_copy=data)

    else:
        # keep OG with paralog
        # choose paralog and get OGs
        pass


if __name__ == '__main__':
    main()
