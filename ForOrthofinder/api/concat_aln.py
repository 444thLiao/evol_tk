import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from Bio import AlignIO, SeqIO
from glob import glob
from os.path import join,exists,dirname
import click
import random
from collections import defaultdict


def remove_identical_seqs(filename, seed=None):
    if seed is not None:
        random.seed(seed)
    groups_dict = defaultdict(list)
    records = list(SeqIO.parse(filename, format='fasta'))
    ori_length = len(records)
    for record in records:
        groups_dict[str(record.seq)].append(record)
    group_info = open(join(os.path.dirname(filename), "identical_groups.list"), 'w')

    count_ = 0
    new_records = []
    for _, records in groups_dict.items():
        if len(records) == 1:
            new_records.append(records[0])
        else:
            new_records.append(random.choice(records))
            group_name = 'group%s' % count_
            for r in records:
                print("%s\t%s\n" % (group_name, r.id),file=group_info)
            count_ += 1
    if len(new_records) == ori_length:
        print("No identical records found")
    else:
        print("Detect identical records, and now dereplicated into single one. From %s to %s" % (ori_length,
                                                                                                 len(new_records)))
    with open(filename, 'w') as f1:
        SeqIO.write(new_records, f1, format='fasta-2line')


def generate_partition_file(outfile, record_pos_info):
    with open(outfile, 'w') as f1:
        for name, start, end,_ in record_pos_info:
            f1.write(f"Protein, {name} = {start}-{end} \n")

def generate_phy_file(outfile, record_pos_info,genome_ids):
    with open(outfile, 'w') as f1:
        for name, start, end, aln_record in record_pos_info:
            total_num = len(genome_ids)
            num_seq = len(aln_record)
            length_this_aln = aln_record.get_alignment_length()
            f1.write(f'{total_num}        {length_this_aln}\n')
            used_ids = []
            for _ in range(num_seq):
                if aln_record[_,:].id in genome_ids:
                    f1.write(f"{aln_record[_,:].id}        {str(aln_record[_,:].seq)}\n")
                    used_ids.append(aln_record[_,:].id)
            for remained_id in set(genome_ids).difference(set(used_ids)):
                f1.write(f"{remained_id}\n{'-'*length_this_aln}\n")
                
@click.command(help="For concating each aln, if it has some missing part of specific genome, it will use gap(-) to fill it")
@click.option("-i", "indir", help="The input directory which contains all separate aln files")
@click.option("-o", "outfile", default=None,help="path of outfile. default id in the `-i` directory and named `concat_aln.aln`")
@click.option("-s", "suffix", default='aln')
@click.option("-gl", "genome_list", default=None, help="it will read 'selected_genomes.txt', please prepare the file, or indicate the alternative name or path.")
@click.option("-rm_I", "remove_identical", is_flag=True, default=False)
@click.option("-seed", "seed", help='random seed when removing the identical sequences')
@click.option("-ct", "concat_type", help='partition or phy or both',default='partition')
def main(indir, outfile,genome_list, remove_identical, seed,concat_type, suffix='aln'):

    if genome_list is None:
        genome_list = join(indir, 'selected_genomes.txt')
    with open(genome_list, 'r') as f1:
        gids = f1.read().split('\n')
    record_pos_info = []
    gid2record = {gid: '' for gid in gids}

    las_pos = 0
    for idx, aln_file in enumerate(glob(join(indir, '*.%s' % suffix))):
        aln_record = AlignIO.read(aln_file, format='fasta')
        length_this_aln = aln_record.get_alignment_length()
        # record the partition
        name = "part%s" % int(idx + 1)
        start, end = las_pos + 1, length_this_aln + las_pos
        las_pos = end
        record_pos_info.append((name, start, end,aln_record ))
        # done record
        for gid in gid2record:
            records = [_ for _ in aln_record if _.id == gid]
            if records:
                gid2record[gid] += str(records[0].seq)
            else:
                gid2record[gid] += '-' * length_this_aln
    if outfile is None:
        outfile = join(indir, 'concat_aln.aln')
        outpartition = join(indir, 'concat_aln.partition')
        outphy = join(indir, 'concat_aln.phy')
    else:
        if not '/' in outfile:
            outfile = './' + outfile
        if not exists(dirname(outfile)):
            os.makedirs(dirname(outfile))
        outfile = outfile
        outpartition = join(dirname(outfile), 'concat_aln.partition')
        outphy = join(dirname(outfile), 'concat_aln.phy')
        
        
    with open(outfile, 'w') as f1:
        for gid, seq in gid2record.items():
            f1.write(f'>{gid}\n')
            f1.write(f'{seq}\n')
    if remove_identical:
        remove_identical_seqs(outfile, seed=seed)
    if concat_type.lower() in ['both','partition']:
        generate_partition_file(outpartition, record_pos_info)
    if concat_type.lower() in ['both','phy']:
        generate_phy_file(outphy, record_pos_info,gids)
    

if __name__ == '__main__':
    main()
