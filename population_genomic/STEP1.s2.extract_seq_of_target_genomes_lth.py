#!/usr/bin/env python

import io, os
from Bio import SeqIO
import click
from tqdm import tqdm


# extract seqs for target and ref genomes
# it's allowed that ref genome is not one of the target genomes


def filter_genomes(records, genome_set):
    new_records = []
    records = sorted(records,key=lambda x:x.description.split('|')[1])
    for record in records:
        header = record.description
        genome = header.split('|')[1]
        if genome in genome_set:
            new_records.append(record)
    return new_records


# xmfa = "mauve_out.xmfa.renamed"
# genome_files = ["target_genomes.txt", "ref_genome.txt"]


@click.command()
@click.option("-i", "xmfa", required=False, default="./mauve_out.xmfa.renamed")
@click.option("-g", "genome_files", multiple=True, required=False, default=["target_genomes.txt", "ref_genome.txt"])
@click.option("-o", "output_f", required=False, default="./mauve_out.xmfa.renamed.selected")
def main(xmfa, output_f, genome_files):
    xmfa = os.path.abspath(xmfa)
    output_f = os.path.abspath(output_f)

    os.makedirs(os.path.dirname(output_f), exist_ok=True)

    list_genomes = set()
    for ffile in genome_files:
        with open(ffile) as f:
            list_genomes.update(set([line[:-1].split('\t')[0]
                                     for line in f
                                     if not line.startswith("#") and len(line) > 1]))
    list_genomes = sorted(list(list_genomes))

    with open(xmfa) as f, open(output_f, "w") as f2:
        segment_list = f.read().split('=\n')
        for segment in tqdm(segment_list):
            records = SeqIO.parse(io.StringIO(segment),
                                  format='fasta')
            new_records = filter_genomes(records, list_genomes)
            SeqIO.write(new_records,
                        handle=f2,
                        format='fasta-2line')
            if segment:
                f2.write('=\n')


if __name__ == '__main__':
    main()
