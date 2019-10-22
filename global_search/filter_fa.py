from Bio import SeqIO
from os.path import *
from tqdm import tqdm
import os


kofam_scan = '/home-user/thliao/software/kofamscan/exec_annotation'

gene_info = {'kegg': {'nxrA': 'K00370',
                      'nxrB': 'K00371',
                      'hao': 'K10535',
                      'amoA': 'K10944',
                      'amoB': 'K10945',
                      'amoC': 'K10946'},
             'TIGFAM': {'nxrA': '',
                        'nxrB': '',
                        'hao': 'TIGR01703',
                        'amoA': 'TIGR03080',
                        'amoB': 'TIGR03079',
                        'amoC': 'TIGR03078'}}


def filter_fa(in_fa, ofile, gene_name):
    
    run(f"{kofam_scan} -o {ofile} --cpu 64 -f mapper-one-line --no-report-unannotated {infa} -p /home-user/thliao/data/kofam/profiles/{ko}.hmm")
    confirmed_id = [_.strip().split('\t')[0]
                    for _ in open(ofile, 'r').read().split('\n') if _]
    confirmed_id = set(confirmed_id)
    remained_records = [_ for _ in tqdm(SeqIO.parse(
        infa, format='fasta')) if _.id in confirmed_id]
    new_fa_file = join(odir, 'unique_seqs_filtered.faa')
    if exists(filter_id_txt):
        ids = [_.strip() for _ in open(filter_id_txt).read().split('\n')]
        remained_records = [_ for _ in remained_records if (
            _.id not in ids) and (_.id.replace('_', ' ') not in ids)]
    with open(new_fa_file, 'w') as f1:
        SeqIO.write(remained_records, f1, format='fasta-2line')
