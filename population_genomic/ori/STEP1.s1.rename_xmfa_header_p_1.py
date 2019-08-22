#!/usr/bin/env python

import os
import sys

# change the header from "> 1:0-0 + genomic_fasta/1A03794.SPAdes.contig1000.fasta" to ">1|1A03794|0..0|+"


xmfa = "parsnp.xmfa"
if len(sys.argv) > 1:
    xmfa = sys.argv[1]

with open(xmfa) as f, open(xmfa+".renamed", "w") as f2:
    # load sequence ids
    dict_genomes = {} # indexed by number and pointed to genome name; 
    line = f.readline()
    line = f.readline()
    line = f.readline()
    while line.startswith("##SequenceIndex"):
        id = line[:-1].split(" ")[1]
        line = f.readline()
        dict_genomes[id] = line.split(" ")[1].split(".")[0]
        f.readline()
        f.readline()
        line = f.readline()
    line = f.readline()
    # load segments one by one
    while line:
        while line and not line.startswith("="):
            id = line[1:].split(":")[0]
            genome = dict_genomes[id]
            pos = line.split(":")[1].split(" ")[0].replace("-", "..")
            strand = line.split(" ")[1]
            header = "%s|%s|%s|%s" % (id, genome, pos, strand)
            seq = ""
            line = f.readline()
            while line and not line.startswith(">") and not line.startswith("="):
                seq += line[:-1]
                line = f.readline()
            f2.write(">%s\n%s\n" % (header, seq))
        f2.write("=\n")
        line = f.readline()

	
    
