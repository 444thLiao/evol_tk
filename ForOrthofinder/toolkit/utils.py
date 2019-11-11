from Bio import SeqIO
import pandas as pd
from subprocess import check_call
import os


def get_dict(file):
    contents = open(file).read().split('\n')
    return_dict = {}
    rev_dict = {}
    for _ in contents:
        id, sep, name = _.partition(":")
        full_name = name.strip()
        if '.' in full_name:
            full_name = full_name.rpartition('.')[0] 
        # extract name precisely, for `6: GCA_000011965.2.faa`, it need to remove suffix and get 'GCA_000011965.2'
        # for `0_9: 000007985v2_00010 Phytochrome-like protein cph1`, it need to get '000007985v2_00010'
        # if not '.' in full_name, it will return ['','',full_name]
        name = full_name.split(' ')[0]
        if name and id:
            return_dict[id] = full_name
            rev_dict[name] = id
    return return_dict, rev_dict


def get_summary_statistic(SC_data):
    number_genomes_presence = SC_data.count(1).sort_values(ascending=False)
    return number_genomes_presence


def get_protein(genomes_path, protein_id):
    if not os.path.exists(genomes_path):
        print(genomes_path, ' not exists')
    all_fa = SeqIO.parse(genomes_path, format='fasta')
    for _ in all_fa:
        if _.id.strip() == protein_id:
            return _
    return None


def get_single_copy(infile):
    data = pd.read_csv(infile, sep='\t', index_col=0, low_memory=False)
    single_copy_mask_df = data.applymap(lambda x: (',' not in str(x))
                                                  and
                                                  (not pd.isna(x)))
    # not contains COMMA(,) and is not NaN
    single_copy_data = data[single_copy_mask_df]
    no_Scopy_genomes = single_copy_mask_df.columns[single_copy_mask_df.isna().all(0)]
    no_Scopy_OG = single_copy_mask_df.index[single_copy_mask_df.isna().all(1)]

    single_copy_data = single_copy_data.loc[single_copy_data.index.difference(no_Scopy_OG),
                                            single_copy_data.columns.difference(no_Scopy_genomes)]
    return single_copy_data


def run_cmd(cmd,**kwargs):
    check_call(cmd,
               shell=True,
               stderr=open('/dev/null','w'),
               stdout=open('/dev/null','w'),
               **kwargs)
