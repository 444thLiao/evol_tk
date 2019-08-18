#!/usr/bin/env python
from __future__ import print_function
from Bio import SeqIO
import os, sys
import click

# validate the phase file by checking the dna extract from the original genome file of the ref genome

def comple(dna):
    dna_comple = {"A": "T",
                  "T": "A",
                  "C": "G",
                  "G": "C"}
    ret = dna_comple.get(dna, "?")
    if ret == "?":
        raise Exception("Error dna: %d" % dna)
    return ret


def validate(filename, contig_name, log_file=sys.stdout):
    f = open(filename, 'r')
    # load seqs from output phase file
    flag = 0
    for nrow, line in enumerate(f):
        line = line[:-1]
        if nrow == 0:
            num_ind = int(line)
            continue
        elif nrow == 1:
            num_SNPs = int(line)
            continue
        elif nrow == 2:
            list_SNPs = line.split(" ")[1:]
            continue
        if list_genomes[nrow - 3] != ref_genome:
            continue
        if num_SNPs != len(list_SNPs):
            mssg = "[Error]: inconsistent num_SNPS=%d and len(list_SNPs)=%d in %s !\n" % (num_SNPs, len(list_SNPs), contig_name)
            print(mssg, file=log_file)

        list_SNPs = [int(list_SNPs[i])
                     for i in range(num_SNPs)]
        # take the number of required num_SNPs as input

        ref_SNPs = line[::]
        # exam the sequence by extracting dna from contig seq of reference genome
        ref_SNPs0 = ""
        ref_SNPs_plus1 = ""
        contig_seq = dict_contigs[contig_name]
        for i in range(num_SNPs):
            pos = list_SNPs[i]
            ref_SNPs0 += str(contig_seq[pos - 1])
            ref_SNPs_plus1 += str(contig_seq[pos])
        if ref_SNPs0 == ref_SNPs:
            print(contig_name, list_genomes[nrow - 3], '1-coordinate position information detected', file=log_file)
            flag = 1
        elif ref_SNPs == ref_SNPs_plus1:
            print(contig_name, list_genomes[nrow - 3], '0-coordinate position information detected', file=log_file)
            flag = 1

    return flag

@click.command()
@click.option("-log","log",default="log.STEP1.s4.txt",required=False)
def main(log):
    log_file = open(log, "w")
    for contig_name in list_contigs:
        ffile = "bialle_SNP.ref_%s.phase" % contig_name
        if ref_not_target:
            ffile = "bialle_SNP.ref_%s.phase.validation" % contig_name
        if not os.path.exists(ffile):
            print(ffile, "not exist", file=log_file)
            continue
        flag = validate(filename=ffile, contig_name=contig_name, log_file=log_file)
        if not flag:
            print("[Error]: Inconsistent SNPs in %s: \n" % (contig_name), file=log_file)
        else:
            print("Validation passed for %s!\n" % contig_name, file=log_file)


if __name__ == '__main__':
    list_genomes = []
    with open("target_genomes.txt") as f:
        for line in f:
            list_genomes.append(line[:-1].split("\t")[0])

    with open("ref_genome.txt") as f:
        ref_genome = f.readline()[:-1]

    records = SeqIO.parse("./%s.fasta" % ref_genome, format='fasta')
    dict_contigs = {"%s_%s" % (ref_genome, _.id.split('_')[1]): _
                    for _ in records}
    list_contigs = sorted(list(dict_contigs.keys()))

    ref_not_target = False
    if ref_genome not in list_genomes:
        ref_not_target = True
        list_genomes.append(ref_genome)

    main()
