"""
For summarizing results output by batch_hmm(some script for annotating genes.)
"""
import os
from collections import defaultdict
from glob import glob
from os.path import join, exists,isdir,isfile
import sys
# sys.path.insert(0,'/'.join(__file__.split('/')[:-4]))
import click
import pandas as pd
from tqdm import tqdm

from dating_workflow.step_script import convert_genome_ID_rev, process_path


def retrieve_info(indir, suffix='.tab',test_or_not=False):
    # deal with --tblout instead of --domtblout
    suffix = suffix.strip('.')
    gid2locus2ko = defaultdict(list)
    if isdir(indir):
        files_list = glob(join(indir, f'*.{suffix}'))
    elif isfile(indir):
        files_list = [indir]
    else:
        raise Exception()
    
    if test_or_not:
        files_list = files_list[:3]
        
    if not files_list:
        exit(f"no files could be found with input {join(indir, f'*.{suffix}')},please check the parameters. ")
    tqdm.write("reading all annotated result")
    for hf in tqdm(files_list):
        for row in open(hf):
            if row.startswith('#'):
                continue
            r = row.split(' ')
            r = [_ for _ in r if _]
            gene_id = r[0]
            ko = r[2]
            evalue = float(r[4])
            gid2locus2ko[convert_genome_ID_rev(gene_id)].append(
                (gene_id, ko, evalue))
    return gid2locus2ko


def filtration_part(gid2locus2ko, evalue=1e-50):
    # filter out with hard threshold of evalue
    post_filtered = {k: [(_[1], _[0], _[2])
                         for _ in v if _[2] <= evalue]
                     for k, v in tqdm(gid2locus2ko.items())}
    # select minimum evalue among all matched KO for each locus
    # TODO: it may be corrected at following version
    ## it could considerate the position overlapping situations
    used_locus = {}
    locus2ko = {}
    tqdm.write("choose best ko for each locus")
    for gid, inlist in tqdm(post_filtered.items()):
        for key, v, evalue in inlist:
            if evalue <= used_locus.get(v, 100):
                used_locus[v] = evalue
                locus2ko[v] = key
    # tqdm.write("choose best ko for each locus")
    post_filtered = defaultdict(lambda: defaultdict(list))
    for locus, ko in locus2ko.items():
        gid = convert_genome_ID_rev(locus)
        post_filtered[gid][ko].append(locus)

    post_filtered = {g: {ko: ','.join(v) for ko, v in d.items()}
                     for g, d in post_filtered.items()}
    return post_filtered


@click.command()
@click.option("-i", "indir", help="input directory.")
@click.option("-o", "odir", help="output directory. If it doesn't exist, it will auto created.")
@click.option("-s", 'suffix', default='tab')
@click.option("-p", 'prefix', default=None, help='prefix of output file, just the file name, it does not need to include dir name. ')
@click.option("-e", "evalue", default=1e-20, help="threshold for filtrations")
@click.option("-t", "transpose", default=False, is_flag=True, help="transpose the output matrix/dataframe or not. default:row is sample/genome, column is KO/annotations")
@click.option("-test", "test", default=False, is_flag=True, help="test the format of the output by sampling only three tabs")
def main(indir, odir, suffix, evalue, transpose, prefix,test):
    indir = process_path(indir)
    odir = process_path(odir)
    gid2locus2ko = retrieve_info(indir, suffix,test_or_not=test)
    post_filtered = filtration_part(gid2locus2ko, evalue)
    if not exists(odir):
        os.makedirs(odir)

    if prefix is not None:
        ofile_info = join(odir, f"{prefix}_info.tab")
        ofile_binary = join(odir, f"{prefix}_binary.tab")
        ofile_num = join(odir, f"{prefix}_num.tab")
    else:
        ofile_info = join(odir, "merged_hmm_info.tab")
        ofile_binary = join(odir, "merged_hmm_binary.tab")
        ofile_num = join(odir, "merged_hmm_num.tab")
    tqdm.write("Complete filterations...")
    tqdm.write("It need time to convert the generated dict into DataFrame. Be patient...")
    final_df = pd.DataFrame.from_dict(post_filtered, orient='index')
    bin_df = final_df.applymap(lambda x: 0 if pd.isna(x) else 1)
    num_df = final_df.applymap(lambda x: 0 if pd.isna(x) else len(str(x).split(',')))
    if transpose:
        final_df = final_df.T
        bin_df = bin_df.T
        num_df = num_df.T
    final_df.to_csv(ofile_info, sep='\t', index=1)
    bin_df.to_csv(ofile_binary, sep='\t', index=1)
    num_df.to_csv(ofile_num, sep='\t', index=1)


if __name__ == "__main__":
    main()
