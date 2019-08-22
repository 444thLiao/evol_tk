#!/usr/bin/env python

import os
import sys

# extract seqs for target and ref genomes
# it's allowed that ref genome is not one of the target genomes


xmfa = "mauve_out.xmfa.renamed"
genome_files = ["target_genomes.txt", "ref_genome.txt"]


list_genomes = []
for ffile in genome_files:
  with open(ffile) as f:
    line = f.readline()
    while line:
	if not line.startswith("#") and len(line)>1:
  	    list_genomes.append(line[:-1].split("\t")[0])
  	    #list_genomes.append(line[:-1])
	line = f.readline()
list_genomes = list(set(list_genomes))


with open(xmfa) as f, open("mauve_out.xmfa.renamed.selected", "w") as f2:
    line = f.readline()
    num_seq = 0
    buff = ""
    while line:
	genome = line.split("|")[1]
	header = line
	seq = f.readline()
	if genome in list_genomes:
	    buff += "%s%s" % (header, seq)
	    num_seq += 1
	line = f.readline() 
	if line.startswith("="):
	    if num_seq == len(list_genomes): # only when all desired genomes are presented
	      	f2.write("%s=\n" % buff)
	    num_seq = 0
	    buff = ""
	    line = f.readline()

