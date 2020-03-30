"""
For summarizing results output by batch_interproscan (no script, run it by batch_any)
like 

python3 ~/bin/batch_run/batch_any.py -i ~/data/plancto/rawdata/genome_protein_files -s faa -o ./interpro_scan -ns '_r' -np 10 -cmd "mkdir {ofile};~/software/interproscan-5.38-76.0/interproscan.sh -i {infile} -d {ofile} -cpu 10 -iprlookup -appl TIGRFAM,CDD,Hamap,Pfam "

"""
import os
from collections import defaultdict
from glob import glob
from os.path import join, exists

import click
import pandas as pd
from tqdm import tqdm

from dating_workflow.step_script import convert_genome_ID_rev, process_path

header = ['Protein Accession',
          'Sequence MD5',
          'Sequence Length',
          'Analysis',
          'Signature Accession',
          'Signature Description',
          'Start location',
          'Stop location',
          'Score',
          'Status',
          'Date',
          'InterPro accession',
          'InterPro description']


def retrieve_info(indir,test=False):
    gid2locus2ko = defaultdict(list)
    exists_db = set()
    files_list = glob(join(indir, '*', f'*.tsv'))

    if not files_list:
        exit(
            f"no files could be found with input {join(indir, '*', f'*.tsv')},please check the parameters. ")
    tqdm.write("reading all annotated result")
    if test:
        files_list = files_list[:10]
    for hf in tqdm(files_list):
        for row in open(hf):
            if not row:
                continue
            r = row.split('\t')
            info_dict = dict(zip(header, r))
            gene_id = info_dict['Protein Accession']
            db = info_dict['Analysis']
            sig_id = info_dict['Signature Accession']

            interpro_id = info_dict.get("InterPro accession", '')
            evalue = float(info_dict['Score'])
            Status = info_dict['Status']
            gid2locus2ko[convert_genome_ID_rev(gene_id)].append(
                (gene_id, db, sig_id, interpro_id, evalue, Status))
            exists_db.add(db)
    return gid2locus2ko, exists_db


def _sub_filter(post_filtered, show_progress=True):
    used_locus = {}
    locus2ko = {}
    tqdm.write("choose best ko for each locus")
    if show_progress:
        iter_o = tqdm(post_filtered.items())
    else:
        iter_o = post_filtered.items()
    for _, inlist in iter_o:
        for gene_id, db, sig_id, interpro_id, _evalue in inlist:
            if _evalue <= used_locus.get(gene_id, 100):
                used_locus[gene_id] = _evalue
                locus2ko[gene_id] = (sig_id, db, interpro_id)
    return locus2ko


def filtration_part(gid2locus2ko, exists_db, evalue=1e-50):
    # filter out with hard threshold of evalue
    tqdm.write(f"filter with evalue {evalue}")
    post_filtered = {k: [(gene_id, db, sig_id, interpro_id, _evalue)
                         for (gene_id, db, sig_id, interpro_id, _evalue, Status) in v
                         if _evalue <= evalue]
                     for k, v in tqdm(gid2locus2ko.items())}
    tqdm.write('separating with db, time consumed....')
    sep_db = {_db: {k: [(gene_id, db, sig_id, interpro_id, _evalue)
                        for (gene_id, db, sig_id, interpro_id, _evalue) in v
                        if db == _db]
                    for k, v in post_filtered.items()}
              for _db in tqdm(exists_db)}

    # select minimum evalue among all matched KO for each locus
    # TODO: it may be corrected at following version
    # it could considerate the position overlapping situations
    locus2ko = _sub_filter(post_filtered)
    sep_l2ko = {_db: _sub_filter(v, show_progress=False)
                for _db, v in sep_db.items()}

    return locus2ko, sep_l2ko


def outut_for(l2ko, odir, name='mixed', transpose=False):
    if not exists(odir):
        os.makedirs(odir)
    tqdm.write('converting into locus2gene side by side table...no progress')
    l2ko_df = pd.DataFrame.from_dict(l2ko).T
    if l2ko_df.shape[1] != 3:
        print(f"it might be something wrong for {name}")
        return
    l2ko_df.columns = ["annotated ID", "database", 'interpro ID']
    l2ko_df.loc[:, 'genome'] = [
        convert_genome_ID_rev(_) for _ in l2ko_df.index]
    l2ko_df.to_csv(join(odir, f"{name}_l2ID.tab"),
                   sep='\t', index=1, index_label='locus')

    tqdm.write(f"start to output {name} locus2gene")
    genome2interpro2locus = defaultdict(lambda: defaultdict(set))
    genome2gene2locus = defaultdict(lambda: defaultdict(set))
    for locus, row in tqdm(l2ko_df.iterrows(), total=l2ko_df.shape[0]):
        genome = row['genome']
        gene = row['annotated ID']
        interpro = row['interpro ID']
        genome2gene2locus[genome][gene].add(locus)
        if interpro:
            genome2interpro2locus[genome][interpro].add(locus)
        
    tqdm.write(f"packing......")
    for _, r in enumerate([genome2gene2locus, genome2interpro2locus]):
        if _ == 1:
            fname = f"{name}_interpro"
        else:
            fname = f"{name}"
        ofile_info = join(odir, f"{fname}_info.tab")
        ofile_binary = join(odir, f"{fname}_binary.tab")
        ofile_num = join(odir, f"{fname}_num.tab")
        final_df = pd.DataFrame.from_dict(r, orient='index')
        bin_df = final_df.applymap(lambda x: 0 if pd.isna(x) else 1)
        num_df = final_df.applymap(
            lambda x: 0 if pd.isna(x) else len(str(x).split(',')))
        if transpose:
            final_df = final_df.T
            bin_df = bin_df.T
            num_df = num_df.T
        final_df.to_csv(ofile_info, sep='\t', index=1, index_label='gene')
        bin_df.to_csv(ofile_binary, sep='\t', index=1, index_label='gene')
        num_df.to_csv(ofile_num, sep='\t', index=1, index_label='gene')


@click.command()
@click.option("-i", "indir", help="input directory.")
@click.option("-o", "odir", help="output directory. If it doesn't exist, it will auto created.")
@click.option("-e", "evalue", default=1e-20, help="threshold for filtrations")
@click.option("-t", "transpose", default=False, is_flag=True, help="transpose the output matrix/dataframe or not. default:row is sample/genome, column is KO/annotations")
@click.option("-test", "test", default=False, is_flag=True, help="test")
def main(indir, odir, evalue, transpose,test):
    indir = process_path(indir)
    odir = process_path(odir)
    gid2locus2ko, exists_db = retrieve_info(indir,test=test)
    locus2ko, sep_l2ko = filtration_part(gid2locus2ko, exists_db, evalue)
    tqdm.write("Complete filterations...")
    if not exists(odir):
        os.makedirs(odir)
    outut_for(locus2ko, odir, name='mixed', transpose=transpose)
    for db, l2ko in sep_l2ko.items():
        outut_for(l2ko,
                  join(odir, f'annotated_with_{db}'),
                  name=db,
                  transpose=transpose)


if __name__ == "__main__":
    main()
