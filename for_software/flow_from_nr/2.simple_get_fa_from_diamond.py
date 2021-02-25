import os
from os.path import *

from Bio import SeqIO
from tqdm import tqdm
#
all_seq = './all_seqs.faa'
# get unique faa from original nr sequence file.

if not exists(all_seq):
    a = open('./nr_retrieve_all.diamond').read()
    all_ids = [_.split('\t')[1] for _ in a.split('\n') if _]
    all_ids = set(all_ids)
    print(len(all_ids))
    with open(all_seq, 'w') as f1:
        nr = SeqIO.parse('/home-user/thliao/data/protein_db/NCBI/nr', format='fasta')
        for record in tqdm(nr):
            if record.id in all_ids:
                SeqIO.write(record, f1, format='fasta-2line')
                f1.flush()

def get_seq_by_annotation(in_file, ofile, keyword):
    all_seq = './all_seqs.faa'
    if isinstance(keyword, str):
        match_ids = [row.split('\t')[0]
                     for row in tqdm(open(in_file).read().split('\n'))
                     if keyword in row]
    elif isinstance(keyword, list):
        match_ids = [row.split('\t')[0]
                     for row in tqdm(open(in_file).read().split('\n'))
                     if any(map(lambda x: x in row, keyword))]
    else:
        raise IOError
    match_ids = set(match_ids)
    records = SeqIO.parse(all_seq, format='fasta')
    remained_record = [_ for _ in tqdm(records) if _.id in match_ids]
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    with open(ofile, 'w') as f1:
        SeqIO.write(remained_record, f1, format='fasta')


if __name__ == '__main__':
    get_seq_by_annotation('nr_parse.output', ofile='./nr_retrieve_nxrB/matched_seq.faa', keyword='nxrB')
    get_seq_by_annotation('nr_parse.output', ofile='./nr_retrieve_nxrA/matched_seq.faa', keyword='nxrA')
    get_seq_by_annotation('nr_parse.output', ofile='./nr_retrieve_amoA/matched_seq.faa', keyword=['amoA', 'pmoA', 'pxmA'])
    get_seq_by_annotation('nr_parse.output', ofile='./nr_retrieve_amoB/matched_seq.faa', keyword=['amoB', 'pmoB', 'pxmB'])

    get_seq_by_annotation('nr_parse.output', ofile='./nr_retrieve_hao/matched_seq.faa', keyword='hao')
    get_seq_by_annotation('nr_parse.output', ofile='./nr_retrieve_amoC/matched_seq.faa', keyword=['amoC', 'pmoC', 'pxmC'])
