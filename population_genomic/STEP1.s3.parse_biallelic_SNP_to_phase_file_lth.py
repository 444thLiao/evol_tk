#!/usr/bin/env python

import os
import sys
import pandas as pd
from Bio import AlignIO
import click
from tqdm import tqdm
import io
from collections import defaultdict
import numpy as np


# gaps will be removed; postions are based on contigs in ref genome; one contig one phase file
# biallelic SNPs only
# alignment on the "-" strand of the reference genome will be converted to reversed complementary strands


# generate a dictionary from concat position to its contig and its position in contig
def return_contig(pos, df, coordinate=1):
    boolean_df = ((df.start <= pos) & (df.end >= pos))

    contig_loc = df.index[boolean_df]
    assert len(contig_loc) == 1
    contig_loc = contig_loc[0]
    contig_start, contig_end = df.loc[boolean_df, :].values[0]
    if coordinate:
        # if pos follow 1-coordinate
        contig_pos = pos - contig_start + 1
        # 1-coordinate
    else:
        # pos follow 0-coordinate
        contig_pos = pos - contig_start + 2
        # 1-coordinate
    return contig_loc, contig_pos


def collect_align(mauve_input):
    # split concated alignment into alignment in each contig
    # 1. collect all ref_genome in each segment
    collect_align_record = []
    with open(mauve_input) as f:
        segment_list = f.read().split('=\n')
        segment_list = [_ for _ in segment_list if _]
        for segment in tqdm(segment_list):
            record = AlignIO.read(io.StringIO(segment),
                                  format='fasta')
            collect_align_record.append(record)
    return collect_align_record


def concat_get_alignment(collect_align_record,
                         ref_genome,
                         df,
                         stat,
                         required_num_alleic=2,
                         no_gap=True):
    # 2. convert the pos in segment to pos in contig
    # 2. and concat all segment into one and record the pos
    # presume the segment won't overlap each others
    contig_pos_list = defaultdict(list)
    contig_seg2align = dict()
    tqdm.write("Collecting all alignment and concat them and also convert the position")
    _count = defaultdict(int)
    for align_record in tqdm(collect_align_record):
        ref_genome_record = [seq
                             for seq in align_record
                             if seq.description.split('|')[1] == ref_genome]
        if not ref_genome_record:
            raise Exception("not existing ref genome record")
        ref_seq = ref_genome_record[0]
        # use the ref_seq as the standard
        header = ref_seq.description
        seqid, genome, pos_str, strand = header.split('|')
        start_concat, end_concat = map(int, pos_str.split('..'))
        if strand == '-':
            start_concat, end_concat = end_concat, start_concat
            # reverse_complement the alignment object
            for _record in align_record:
                _record.seq = _record.seq.reverse_complement()

        s_contig, s_pos = return_contig(start_concat, df, coordinate=1)
        e_contig, e_pos = return_contig(end_concat, df, coordinate=1)

        # remove '-' at align_record
        if s_contig != e_contig:
            # this segment spans two or more contigs
            # fixme: not sure how to deal with it
            # keep biggest one? or keep two of them
            print("[Warning]: Start and end of the alignment are from two different contigs! %s and %s " % (s_contig, e_contig))
            #
            pass
        else:

            if contig_seg2align.get(s_contig) is None:
                contig_seg2align[s_contig] = align_record
            else:
                contig_seg2align[s_contig] += align_record

            # init the position number of the alignment
            pos_list = []
            pos = s_pos  # start from s_pos
            for dna in str(ref_seq.seq):
                if dna != '-':
                    pos_list.append(pos)
                    pos += 1
                else:
                    pos_list.append('?')
            assert pos - 1 == e_pos  # assert to make sure no errors
            contig_pos_list[s_contig] += pos_list
        _count[s_contig] += 1
    contig_name = list(contig_pos_list)[0]
    assert len(contig_pos_list[contig_name]) == contig_seg2align[contig_name].get_alignment_length()
    # 3. find biallelic SNPs
    contig2biallelic_SNP_pos = defaultdict(list)
    contig2subset_align = defaultdict(list)
    tqdm.write("start find biallelic SNP")
    for contig_name in tqdm(contig_seg2align):
        align_record = contig_seg2align[contig_name]
        np_align = np.array([list(rec) for rec in align_record])
        num_alleic = np.apply_along_axis(lambda x: len(set(x)), 0, np_align)
        # count the number of alleic, including the gap
        gap_alleic = np.apply_along_axis(lambda x: '-' in x, 0, np_align)
        # True or not contain the gap within one column

        # num_alleic == 2 and gap_alleic == False
        if no_gap:
            condition = (num_alleic == required_num_alleic) & (~gap_alleic)
        else:
            condition = (num_alleic == required_num_alleic)
        _dict = stat[contig_name]
        _dict["num_sites_aligned(w_o_gap)"] = (~gap_alleic).sum()
        _dict["num/contig_size"] = "%.1f%%" % ((100 * (~gap_alleic).sum()) / (_dict["contig_size"]))
        _dict["num_SNPs"] = ((num_alleic >= 2) & (~gap_alleic)).sum()
        _dict["num/aligned"] = "%.1f%%" % (100.0 * _dict["num_SNPs"] / _dict["num_sites_aligned(w_o_gap)"])
        _dict["num_biallelic_SNPs"] = ((num_alleic == 2) & (~gap_alleic)).sum()
        _dict["num/SNPs"] = "%.1f%%" % (100.0 * _dict["num_biallelic_SNPs"] / _dict["num_SNPs"])
        _dict["num_alignment_blocks"] = _count[contig_name]
        _dict["num_individuals"] = np_align.shape[0]

        ori_pos_list = np.array(contig_pos_list[contig_name])
        contig2biallelic_SNP_pos[contig_name] = list(ori_pos_list[condition])

        align_matrix = np_align[:, condition]
        align_text = '\n'.join([''.join(row) for row in align_matrix])
        contig2subset_align[contig_name] = align_text
    return contig2biallelic_SNP_pos, contig2subset_align, stat


# write phase file for each contig
def writePhaseFile(ext, contig2biallelic_SNP_pos, contig2subset_align, ref_genome):
    for name, pos_list in contig2biallelic_SNP_pos.items():
        if not pos_list:
            continue
        with open("bialle_SNP.ref_%s.%s" % (name, ext), "w") as f:
            num_individuals = contig2subset_align[name].count('\n') + 1
            if "validation" in ext:
                num_individuals += 1
            num_SNPs = len(pos_list)
            str_SNP_pos = " ".join(map(str, pos_list))
            f.write("%d\n%d\nP %s\n" % (num_individuals, num_SNPs, str_SNP_pos))
            f.write(contig2subset_align[name])
            if "validation" in ext:
                f.write(contig2subset_align[ref_genome] + "\n")


@click.command()
@click.option("-i", "input_table", default="ref_genome.pos.table", required=False)
@click.option("-i2", "mauve_input", default="mauve_out.xmfa.renamed.selected", required=False)
def main(input_table, mauve_input):
    # prepare a dataframe for collecting the stats
    stat = defaultdict(dict)
    # load contig start and end to Contig object, for ref_genome only
    tqdm.write("start loading required data")
    df = pd.read_csv(input_table,
                     sep='\t',
                     header=None,
                     index_col=0,
                     comment="#")
    df.columns = ["start", 'end']
    for contig_name, row in df.iterrows():
        stat[contig_name]["contig_id"] = contig_name
        stat[contig_name]["contig_size"] = row["end"] - row["start"] + 1
    ref_genome = open(input_table, 'r').readline().strip('\n').split('=')[1]
    collect_align_record = collect_align(mauve_input)

    contig2biallelic_SNP_pos, contig2subset_align, stat = concat_get_alignment(collect_align_record,
                                                                               ref_genome,
                                                                               df,
                                                                               stat)
    stat_df = pd.DataFrame.from_dict(stat, orient='index')
    stat_df = stat_df.reindex(sorted(stat_df.index, key=lambda x: int(x.split('_')[1])))
    stat_df = stat_df.fillna("--")
    with open("log.STEP1.s3.SNP_summary.txt", "w") as f:
        stat_df.to_csv(stat_df, f, index=False)
    writePhaseFile("phase", contig2biallelic_SNP_pos, contig2subset_align, ref_genome)


if __name__ == '__main__':
    main()
