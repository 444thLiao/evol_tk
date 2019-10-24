from Bio import SeqIO
import numpy as np


def get_len(records):
    return [len([bp for bp in _ if bp!='-'])
            for _ in records]


def compared(infile):
    trimal = infile.replace('.aln','.trimal')
    records_ori = list(SeqIO.parse(infile,format='fasta'))
    records_aft = list(SeqIO.parse(trimal,format='fasta'))
    
    ori_medi = np.median(get_len(records_ori))
    ori_mean = np.mean(get_len(records_ori))
    ori_std = np.std(get_len(records_ori))

    aft_medi = np.median(get_len(records_aft))
    aft_mean = np.mean(get_len(records_aft))
    aft_std = np.std(get_len(records_aft))
    
    print(f'ori median: {ori_medi}')
    print(f'ori mean: {ori_mean}')
    print(f'ori std: {ori_std}')
    print(f'aft median: {aft_medi}')
    print(f'aft mean: {aft_mean}')
    print(f'aft std: {aft_std}')
    
    
    
    
    