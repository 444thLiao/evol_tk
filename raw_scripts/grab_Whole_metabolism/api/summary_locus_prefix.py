from glob import glob
import click
from os.path import join, dirname, basename, abspath
import os
from tqdm import tqdm
from Bio import SeqIO
import pandas as pd
import string, random
import gzip

def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_letters
    return ''.join(random.choice(letters) for _ in range(stringLength))


def prokka_summary(indir):
    p_files = glob(join(indir, '*', "*.faa"))
    if not p_files:
        raise Exception("No faa file in %s" % indir)
    sname2locus = {}
    occured_id = []
    for p_file in tqdm(p_files):
        # iterating faa files
        sname = basename(dirname(p_file)).replace('.faa', '')
        records = SeqIO.parse(p_file, format='fasta')
        record = next(records)
        # get the first one, enough
        locus_tag = record.id.split('_')[0]
        if locus_tag in occured_id:
            # if locus_tag occurred multiple times, regenerated a new one
            records = SeqIO.parse(p_file, format='fasta')
            locus_tag = randomString(len(occured_id[-1]) + 1)
            new_records = []
            for record in records:
                record.id = locus_tag + '_' + record.id.split('_')[-1]
                new_records.append(record)
            os.renames(p_file, p_file + '.backup')
            with open(p_file, 'w') as f1:
                f1.write(new_records)
        occured_id.append(locus_tag)
        sname2locus[sname] = {}
        sname2locus[sname]['locus_prefix'] = locus_tag
    result_df = pd.DataFrame.from_dict(sname2locus, orient='index')
    return result_df


def download_summary(indir):
    genome_dirs = glob(join(indir, '*', '*'))
    genome_dirs = [_ for _ in genome_dirs
                   if os.path.isdir(_)]
    if not genome_dirs:
        raise Exception("No directory detected in %s" % indir)
    missing_faa_samples = []
    sname2locus = {}
    occured_id = []
    for each_dir in tqdm(genome_dirs):
        # iterating all directory
        p_faa = glob(join(each_dir, '*.faa.gz'))
        sname = basename(each_dir)
        if not p_faa:
            # if not download faa file, pass it and record it.
            missing_faa_samples.append(sname)
            sname2locus[sname] = {}
            continue
        p_faa = p_faa[0]
        records = SeqIO.parse(gzip.open(p_faa,'rt'), format='fasta')
        random_prefix = randomString(10)
        if random_prefix in occured_id:
            random_prefix = randomString(10)
        occured_id.append(random_prefix)
        new_records = []
        for record in records:
            record.id = random_prefix + '_' + record.id
            new_records.append(record)
        with open(join(each_dir, 'generated_protein.faa'), 'w') as f1:
            SeqIO.write(new_records, f1, format='fasta-2line')
        sname2locus[sname] = {}
        sname2locus[sname]['locus_prefix'] = random_prefix
    result_df = pd.DataFrame.from_dict(sname2locus, orient='index')
    return result_df


@click.command(help="quickly get a summary file from prokka_o")
@click.option("-i", "indir", help='input dir, normally is the output directory of prokka.')
@click.option("-o", "outfile", help='output summary file')
@click.option("-t", "typeOfdata", help='data type including prokka or download')
def main(indir, outfile, typeOfdata):
    indir = abspath(indir)

    if not os.path.exists(dirname(outfile)):
        os.makedirs(dirname(outfile), exist_ok=True)
    if typeOfdata.lower() == 'prokka':
        result_df = prokka_summary(indir)
    elif typeOfdata.lower() == 'download':
        result_df = download_summary(indir)
    else:
        raise Exception('accepted parameters of -t included')
    result_df.to_csv(outfile, index_label='sample_name')

if __name__ == '__main__':
    main()