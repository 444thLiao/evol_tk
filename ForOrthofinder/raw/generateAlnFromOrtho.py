import pandas as pd
from Bio import SeqIO
from os.path import join, dirname, abspath

indir = './data/Results_Aug09_1'


orthogroups_path = join(indir, "Orthogroups.csv")


def get_single_copy(infile):
    data = pd.read_csv(infile, sep='\t', index_col=0)
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


def get_summary_statistic(SC_data):
    number_genomes_presence = SC_data.count(1).sort_values(ascending=False)
    return number_genomes_presence


def get_clustered_gene_OG(df, column):
    sub_df = df.loc[:, column]
    idx = sub_df.index[~sub_df.isna()]
    if len(idx) != 1:
        print("Multiple OG detected? Please verify the column provided it correct")
    if len(idx) != 0:
        return idx[0]

def get_clustered_genomes(df, column):
    sub_df = df.loc[:, column]
    idx = df.index[~sub_df.isna()]
    sub_df = df.loc[idx[0], :]
    genomes = df.columns[~sub_df.isna()]
    if len(genomes) != 0:
        return genomes


def get_dict(file):
    contents = open(file).read().split('\n')
    return_dict = {}
    rev_dict = {}
    for _ in contents:
        id, sep, name = _.partition(":")
        full_name = name.strip()
        name = full_name.split('.')[0]
        name = name.split(' ')[0]
        if name and id:
            return_dict[id] = full_name
            rev_dict[name] = id
    return return_dict, rev_dict


def get_protein(genomes_path, protein_id):
    all_fa = SeqIO.parse(genomes_path, format='fasta')
    for _ in all_fa:
        if _.id.strip() == protein_id:
            return _
    return None


def get_seq_with_OG(orthogroups_path, OG, output_dir, single_copy=True):
    if type(OG) == str:
        OG = [OG]
    thisdir = dirname(orthogroups_path)
    SeqID_file = join(thisdir, "SequenceIDs.txt")
    SpeID_file = join(thisdir, "SpeciesIDs.txt")
    id2seq, seq2id = get_dict(SeqID_file)
    id2spe, spe2id = get_dict(SpeID_file)
    if single_copy:
        data = get_single_copy(orthogroups_path)
    else:
        data = pd.read_csv(orthogroups_path, sep='\t', index_col=0)
    if set(OG).difference(set(data.index)):
        raise Exception("Some OG is not presented at the index of data")
    sub_data = data.loc[OG, :]
    sub_data = sub_data.loc[:, ~sub_data.isna().all(0)]
    # remove all nan genomes
    # open these genomes
    sub_data.columns = [spe2id[_] for _ in sub_data.columns]
    sub_data = sub_data.applymap(lambda x: seq2id.get(x,'nan'))
    species_path_temp = join(thisdir, "Species{speid}.fa")

    for og in OG:
        with open(join(output_dir, og + '.faa'), 'w') as f1:
            seqs = []
            for seq_id in sub_data.loc[og, :]:
                if seq_id == 'nan':
                    continue
                speid = seq_id.split('_')[0]
                spe_file = species_path_temp.format(speid=speid)
                record = get_protein(spe_file, seq_id)
                if record is not None:
                    record.id = id2seq[record.id]
                    seqs.append(record)
                # else:
                    # print(seq_id, "doesn't exist?")
            SeqIO.write(seqs, handle=f1, format='fasta-2line')


if __name__ == '__main__':
    single_copy_df = get_single_copy(orthogroups_path)
    number_genomes_presence = get_summary_statistic(single_copy_df)
