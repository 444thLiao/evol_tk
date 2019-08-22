import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from Bio import AlignIO
from glob import glob
from os.path import join
import click


def generate_partition_file(outfile, record_pos_info):
    with open(outfile, 'w') as f1:
        for name, start, end in record_pos_info:
            f1.write(f"Protein, {name} = {start}-{end} \n")


@click.command(help="For concating each aln, if it has some missing part of specific genome, it will use gap(-) to fill it")
@click.option("-i", "indir",help="The input directory which contains all separate aln files")
@click.option("-s", "suffix", default='aln')
@click.option("-gl", "genome_list", default=None,help="it will read 'selected_genomes.txt', please prepare the file, or indicate the alternative name or path.")
def main(indir, genome_list, suffix='aln'):
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
        record_pos_info.append((name, start, end))
        # done record
        for gid in gid2record:
            records = [_ for _ in aln_record if _.id == gid]
            if records:
                gid2record[gid] += str(records[0].seq)
            else:
                gid2record[gid] += '-' * length_this_aln
    with open(join(indir, 'concat_aln.aln'), 'w') as f1:
        for gid, seq in gid2record.items():
            f1.write(f'>{gid}\n')
            f1.write(f'{seq}\n')
    generate_partition_file(join(indir, 'concat_aln.partition'), record_pos_info)


if __name__ == '__main__':
    main()
