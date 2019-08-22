#!/usr/bin/env python

import os
import sys

# validate the phase file by checking the dna extract from the original genome file of the ref genome


# list of genome names in the same order as phase file
list_genomes = []
with open("target_genomes.txt")  as f:
    line = f.readline()
    while line:
	list_genomes.append(line[:-1].split("\t")[0])
	line = f.readline()


# get sequence of ref genome
dict_contigs = {} # indexed by contig name and pointed to contig seq
ref_genome = "" 
list_contigs = []
with open("ref_genome.txt") as f:
    ref_genome = f.readline()[:-1]
with open("./%s.fasta" % ref_genome) as f:
    line = f.readline()
    while line:
	contig_name = "%s_%s" % (ref_genome, line.split("_")[1].split(" ")[0])
	dict_contigs[contig_name] = f.readline()[:-1]
	list_contigs.append(contig_name)
	line = f.readline()
list_contigs = sorted(list_contigs)

ref_not_target = 0
if ref_genome not in list_genomes:
    ref_not_target = 1
    list_genomes.append(ref_genome)

def comple(dna):
    ret = "?"
    if dna=="A":
	ret ="T"
    elif dna=="T":
	ret = "A"
    elif dna=="G":
	ret = "C"
    elif dna =="C":
	ret = "G"
    else:
	print "Error dna: %d" % dna
	sys.exit()
    return ret

# validate phase file content
summ = ""
for contig_name in list_contigs:
    ffile = "bialle_SNP.ref_%s.phase" % contig_name
    if ref_not_target:
	ffile = "bialle_SNP.ref_%s.phase.validation" % contig_name
    print ffile
    if os.path.isfile(ffile):
      with open(ffile) as f:
	# load seqs from output phase file
	num_ind = int(f.readline()[:-1])
	num_SNPs = int(f.readline()[:-1])
	#print num_SNPS
	list_SNPs = f.readline()[:-1].split(" ")[1:]
	if num_SNPs != len(list_SNPs):
	    mssg = "[Error]: inconsistent num_SNPS=%d and len(list_SNPs)=%d in %s !\n" % (num_SNPs, len(list_SNPs), contig_name)
	    summ += mssg
	for i in range(num_SNPs):
	    list_SNPs[i] = int(list_SNPs[i])
	ref_SNPs = "" # list of SNPs of reference genome
	i_genome = 0
	line = f.readline()
	while line:
	    if list_genomes[i_genome] == ref_genome:
	        ref_SNPs = line[:-1]
		break
	    i_genome += 1
	    line = f.readline()

	# exam the sequence by extracting dna from contig seq of reference genome
	ref_SNPs0 = ""
	ref_SNPs_plus1 = ""
	ref_SNPs_minus1 = ""
	contig_seq = dict_contigs[contig_name]
	for i in range(num_SNPs):
	    pos = list_SNPs[i]
	    ref_SNPs0 += contig_seq[pos-1]
	    ref_SNPs_plus1  += contig_seq[pos]
	    ref_SNPs_minus1 += contig_seq[pos-2]
	flag = 0
	for i in range(num_SNPs):
	    #if ref_SNPs[i] != ref_SNPs0[i] and ref_SNPs[i]!=comple(ref_SNPs0[i]):
	    if ref_SNPs[i] != ref_SNPs_plus1[i]:# and ref_SNPs[i]!=comple(ref_SNPs_plus1[i]):
		flag = 1
		break
	if flag:
	    summ +=  "[Error]: Inconsistent SNPs in %s: \nphase file = \n'%s...'\nextracted from genome = \n'%s...'\n" %(contig_name, ref_SNPs[:1000], ref_SNPs0[:1000]) 
	    summ +=  "'%s'\n'%s\n'" % (ref_SNPs_plus1, ref_SNPs_minus1)
	else:
	    summ +=  "Validation passed for %s!\n" % contig_name 
	summ += "\n"

print summ
with open("log.STEP1.s4.txt", "w") as f:
    f.write(summ)
