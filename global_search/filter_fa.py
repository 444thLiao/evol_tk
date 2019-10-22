from Bio import SeqIO
from os.path import *
from tqdm import tqdm
import os
import numpy as np
from subprocess import check_call
import sys
def run(cmd):
    check_call(cmd,shell=True,stdout=open('/dev/null'))
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


def filter_fa_by_db(infa, ofile, gene_name):
    _ginfo = gene_info['kegg']
    ko = _ginfo.get(gene_name,'')
    if not ko:
        return
    otab=join(dirname(infa),f'{ko}.kofam.out')
    if not exists(otab):
        run(f"{kofam_scan} -o {otab} --cpu 64 -f mapper-one-line --no-report-unannotated {infa}")
    confirmed_id = [_.strip().split('\t')[0]
                    for _ in open(otab, 'r').read().split('\n') 
                    if _ and ko in _]
    confirmed_id = set(confirmed_id)
    remained_records = [_ 
                        for _ in tqdm(SeqIO.parse(infa, format='fasta')) 
                        if _.id in confirmed_id]
    with open(ofile, 'w') as f1:
        SeqIO.write(remained_records, f1, format='fasta-2line')
    return ofile

def filter_fa_by_length_dis(in_fa, ofile=None,output_records=True,down_threshold=25,upper_threshold=100,hard_filter=None):
    records = [_ for _ in SeqIO.parse(in_fa,format='fasta')]
    length_dis = [len(_.seq) for _ in records]
    down_len = np.percentile(length_dis,down_threshold)
    upper_len = np.percentile(length_dis,upper_threshold)
    print('down len: ',down_len)
    print('upper len: ',upper_len)
    _s = sorted(records,key=lambda x:len(x.seq))
    print('longest seq is %s, has %s AA' % (_s[-1].id,len(_s[-1].seq)))
    print('shortest seq is %s, has %s AA' % (_s[0].id,len(_s[0].seq)))
    if hard_filter is None:
        remained_records = [_
                        for _ in records
                        if len(_.seq)>down_len and len(_.seq)<upper_len]
    else:
        remained_records = [_
                        for _ in records
                        if len(_.seq)>hard_filter]
    print('ori number of sequences: ',len(records))
    print('remained number of sequences: ',len(remained_records))
    _s = sorted(remained_records,key=lambda x:len(x.seq))
    print('longest seq is %s, has %s AA' % (_s[-1].id,len(_s[-1].seq)))
    print('shortest seq is %s, has %s AA' % (_s[0].id,len(_s[0].seq)))
    if (ofile is not None) and (not output_records):
        if not exists(dirname(ofile)):
            os.makedirs(dirname(ofile))
        with open(ofile,'w') as f1:
            
            SeqIO.write(remained_records,f1,format='fasta-2line')
    elif output_records:
        return remained_records
    return 

def cluster_fa(infa):
    if not '/' in infa:
        infa = './' + infa
    odir = dirname(infa)
    for threshold in [90,95,98]:
        run(f"cd-hit -i {infa} -o {odir}/cluster_{threshold} -c 0.{threshold} -n 5  -T 20 -d 0")
        
if __name__ == "__main__":
    for gene in ['nxrA', 'nxrB', 'amoA', 'amoB']:
        print('processing !! ',gene)
        odir = f'./nr_retrieve_{gene}'
        infa = join(odir,'matched_seq.faa')
        new_fa = filter_fa_by_db(infa,
                                         join(odir,'filtered_by_kegg.faa'),
                                         gene)
        records = [_ for _ in SeqIO.parse(new_fa,format='fasta')]
        if len(records) > 1500:
            cluster_fa(infa)
            final_clustered_fa = join(odir,'cluster_98')
            final_fa = final_clustered_fa +'_filtered_lengths.fa'
            # final_fa = join(odir,'filtered_by_length.faa')
            
            if gene='amoB':
                hard_f_value = 240
            elif gene = 'amoA':
                final_clustered_fa = join(odir,'cluster_95')
                final_fa = final_clustered_fa +'_filtered_lengths.fa'
                hard_f_value = 200
            filter_fa_by_length_dis(final_clustered_fa,
                                    final_fa,
                                    output_records=False,
                                    hard_filter=hard_f_value,)
                                    #down_threshold=down_threshold)
            # records = [_ for _ in SeqIO.parse(final_fa,format='fasta')]
            # if len(records) >= 1000:
            #     cluster_fa(final_fa)
        
