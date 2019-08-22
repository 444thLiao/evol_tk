#!/usr/bin/env python
from __future__ import print_function
import os
import io
import click
from Bio import SeqIO
from tqdm import tqdm


# change the header from "> 1:0-0 + genomic_fasta/1A03794.SPAdes.contig1000.fasta" to ">1|1A03794|0..0|+"

def collect_comment(comment):
    seq_indexs = [_.split(" ")[1]
                  for _ in comment
                  if _.startswith("##SequenceIndex")]
    seq_headers = [_.split(" ")[1].split(".")[0]
                   for _ in comment
                   if _.startswith("##SequenceFile")]
    dict_genomes = dict(zip(seq_indexs,
                            seq_headers))
    return dict_genomes


def rename(records, id2genome_name):
    new_records = []
    for record in records:
        before_header = record.description
        id = record.id.split(":")[0]
        genome = id2genome_name[id]
        pos = record.id.split(":")[1].split(" ")[0].replace("-", "..")
        strand = before_header.split(' ')[1]
        new_header = "%s|%s|%s|%s" % (id,
                                      genome,
                                      pos,
                                      strand)
        record.id = record.name = record.description = ''
        record.id = new_header
        record.seq = record.seq.upper()
        new_records.append(record)
    return new_records


@click.command()
@click.option("-i", "xmfa", default="./parsnp.xmfa", required=False)
def main(xmfa):
    xmfa = os.path.abspath(xmfa)
    input_f = open(xmfa, 'r')
    output_f = open(os.path.join(os.path.dirname(xmfa),
                                 "mauve_out.xmfa.renamed"),
                    "w")

    # load sequence ids
    input_contents = input_f.read().split("\n")
    comment_contents = [_ for _ in input_contents if _.startswith("#")]
    seqs_contents = [_ for _ in input_contents if not _.startswith("#")]

    id2genome_name = collect_comment(comment_contents)  # indexed by number and pointed to genome name;

    segment_list = [segment
                    for segment in '\n'.join(seqs_contents).split('=\n')]

    for segment in tqdm(segment_list):
        records = SeqIO.parse(io.StringIO(segment), format='fasta')
        new_records = rename(records, id2genome_name)
        SeqIO.write(new_records,
                    handle=output_f,
                    format='fasta-2line')
        if segment:
            output_f.write("=\n")
    output_f.close()
    input_f.close()


if __name__ == '__main__':
    main()
