#!/usr/bin/env python

import os
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# gaps will be removed; postions are based on contigs in ref genome; one contig one phase file
# biallelic SNPs only
# alignment on the "-" strand of the reference genome will be converted to reversed complementary strands


# Contig of reference genome
class Contig():

    def __init__(self, arr):
 	self.start = int(arr[1])	# start of contig as if all contigs are concatenated
	self.end = int(arr[2])		# end of contig as if all contigs are concatenated
	self.name = arr[0]		# contig name, e.g. 1A01234_2
	# orignal alignment
	self.dict_alignments = {} 	# dict of mauve alignments, each of which is indexed by the alignment start and pointed to  
				      		# a dict_seq who is indexed by genome name and pointed to the sequence
				      		# e.g. dict_alignment[211] = dict_seq; dict_seq["1A01234"]= "ATCCGG..."
	self.list_alignment_starts = [] # list of alignment starts that are corresponding to dict_alignments
	# SNP
	self.list_SNP_pos = []	# list of biallelic SNP positions
	self.dict_SNPs = {}	# dictionary of biallelic SNPs, indexed by genome name and pointed to the SNP corresponding to list_SNP_pos
	self.num_multi_locus_SNPs = 0 # num of multi-locus SNP positions (>2 alleles)i, for stat only
	self.num_aligned_sites = 0 	# number of aligned sites that don't have any gap

    # convert position as if all contigs are concatenated to the real position in the contig
    def getPosInContig(self, pos):
	ret = -1
	if pos >= self.start and pos <= self.end:
	    ret = pos - self.start + 1
	return ret



# load contig start and end to Contig object, for ref_genome only
dict_contigs = {} 	# dict of Contig objects
list_contig_inds = [] 	# list of contig names, suppose the contig name is in form of e.g. "1A01234_12"
ref_genome = ""		
with open("ref_genome.pos.table") as f:
    line = f.readline()
    ref_genome = line[:-1].split("=")[1]
    while line:
	if not line.startswith("#"):
	    arr = line[:-1].split("\t")
	    contig = Contig(arr)
	    dict_contigs[contig.name] = contig
	    list_contig_inds.append(int(contig.name.split("_")[1]))
	line = f.readline()
list_contig_inds = sorted(list_contig_inds)
list_contig_names = []
for i in range(len(list_contig_inds)):
    list_contig_names.append("%s_%d" % (ref_genome, list_contig_inds[i]))


# load contig.list_alignments for each contig in the ref genome
with open("mauve_out.xmfa.renamed.selected") as f:
    line = f.readline()
    while line:
	# 1) load alignment sequence
	dict_seq = {} # indexed by genome name and pointed to the seq in the alignment
	ref_contig = ""   # the alignment belongs to which contig
	start_in_contig = 0  # the start position of the alignment in the contig
	end_in_contig = 0    # the end position of the alignment in the contig
	strand_in_contig = "+"	# which strand of the contig the alignment is accouted for
	adjust_align_start = 0  # num of bp cut from the start of the alignment (NOT CONTIG!!!); the alignment will be cut to [adjust_align_start:adjust_align_end]
	adjust_align_end = 0	# num of bp cut form the end of the alignment (NOT CONTIG!!!) (<0)
	while line and not line.startswith("="):
	    genome = line.split("|")[1]#.split("_")[2]
	    pos_str = line.split("|")[2]
	    strand = line.split("|")[3][0]
	    seq = f.readline()[:-1]
	    dict_seq[genome] = seq
	    # record pos_in_ref_genome if ref_genome
	    if genome == ref_genome:
	      strand_in_contig = strand
              start_in_cat = int(pos_str.split(".")[0]) # pos in concatinated seq
              end_in_cat = int(pos_str.split(".")[2])
	      # if in the "+" strand
	      if strand == "+":
		# find position in contig
		for i in range(len(list_contig_names)):
		    name = list_contig_names[i]
		    contig = dict_contigs[name]
		    start = contig.getPosInContig(start_in_cat) # pos in contig
		    end = contig.getPosInContig(end_in_cat)
		    if start > 0:
			start_in_contig = start
			ref_contig = name
			# the end position of the alignment in the current contig is found
			if end > 0:
			    end_in_contig = end
			    break
			# else, the alignment is involved in two contigs: must select the contig where the major part of alignment is located.
			else:
			    next_contig = dict_contigs[list_contig_names[i+1]]
			    print "[Warning]: Start and end of the alignment are from two different contigs! %s and %s " % (contig.name, next_contig.name)
			    print "           contig_%s_in_cat=%d..%d   aln_in_cat=%s   aln_distr_to_contig_end=[%d,%d]   alignment_stradn=%s "\
				 % (name, contig.start, contig.end, pos_str, contig.end-start_in_cat, contig.end-end_in_cat, strand_in_contig)
			    seq_len = len(dict_seq[ref_genome].replace("-", ""))
			    # the major part is in current contig, then cut the alignment
			    if contig.end-start_in_cat >= 0.5*seq_len:
				adjust_align_end = (contig.end - end_in_cat)
				end_in_contig = contig.end-contig.start+1 # end of the contig
				break
			    # the major part is in another contig, go to the next contig
                            else:
                                adjust_align_start = (contig.end-start_in_cat) + 1
                                start_in_cat = contig.end+1
                                continue
	      # if on the "-" strand
	      else:
		# find position in contig
                for i in range(len(list_contig_names)):
                    name = list_contig_names[len(list_contig_names) - i -1] # scan the contig from the opposite direction
                    contig = dict_contigs[name]
                    start = contig.getPosInContig(start_in_cat) # pos in contig
                    end = contig.getPosInContig(end_in_cat)
                    if end > 0:		# find end first, because on the opposite strand
                        end_in_contig = end
			ref_contig = name
                        # the start position of the alignment in the current contig is found
                        if start > 0:
                            start_in_contig = start
                            break
                        # else, the alignment is involved in two contigs: must select the contig where the major part of alignment is located.
                        else:
			    next_contig = dict_contigs[list_contig_names[len(list_contig_names)-i-2]]
                            print "[Warning]: Start and end of the alignment are from two different contigs! %s and %s" % (contig.name, next_contig.name)
                            print "           contig_%s_in_cat=%d..%d   aln_in_cat=%s   aln_distr_to_contig_start=[%d,%d]   alignment_strand=%s"\
                                 % (name, contig.start, contig.end, pos_str, start_in_cat-contig.start, end_in_cat-contig.start, strand_in_contig)
                            seq_len = len(dict_seq[ref_genome].replace("-", ""))
                            # the major part is in current contig, then cut the alignment
                            if end_in_cat-contig.start+1 >= 0.5*seq_len:
                                adjust_align_start = (contig.start - start_in_cat)
				start_in_contig = 1
                                break
                            # the major part is in another contig, go to the next contig
                            else:
                                adjust_align_end = (contig.start-end_in_cat)-1
                                end_in_cat = contig.start-1
                                continue
	    line = f.readline()
	# if the alignment is on the "-" strand of the ref genome, reverse complement the whole alignment
   	if strand_in_contig == "-":
	    for g,seq in dict_seq.items():
		dna = Seq(seq, generic_dna)
		dict_seq[g] = str(dna.reverse_complement())
	# adjust alignment start/end based on adjust_seq
	if adjust_align_start > 0:
	    # find the length of adjust_align_start without gap in ref_genome
	    pos = 0
	    seq = dict_seq[ref_genome]
	    num = 0
	    while num <= adjust_align_start:
		if seq[pos] != "-":
		    num+= 1
		pos += 1
	    pos -= 1
	    for g,seq in dict_seq.items():
		dict_seq[g] = seq[pos:]
	    print "      Start of alignment is vertically cut from %d.\n" % (pos)
	if adjust_align_end < 0:
	    # find the length of adjust_align_end without gap in ref_genome
	    pos = 0
	    num = 0
	    seq = dict_seq[ref_genome]
	    while num <= abs(adjust_align_end):
		if seq[pos]!="-":
		    num+=1
		pos -= 1
	    for g,seq in dict_seq.items():
                dict_seq[g] = seq[:(pos+1)]
	    print "      End of alignment is vertically cut to %d\n" % (pos)
	# make a separate phase file for one contig
	if start_in_contig>0 and end_in_contig>0:
	    contig = dict_contigs[ref_contig]
	    # record alignment
	    contig.dict_alignments[start_in_contig] = dict_seq
	    contig.list_alignment_starts.append(start_in_contig)	
	    # initiate SNP dict
	    for g in dict_seq:
		contig.dict_SNPs[g] = ""
	# prepare for next alignment
	line = f.readline()
	
# find biallelic SNPs
for name in list_contig_names:
    # for each contig in ref genome
    contig = dict_contigs[name]
    for i in range(len(contig.list_alignment_starts)):
	# for each alignment in the contig
	align_start = contig.list_alignment_starts[i]
	dict_seq = contig.dict_alignments[align_start]
	ref_seq = dict_seq[ref_genome]
	for j in range(len(dict_seq[ref_genome])):
	    # check biallelic snp for each site
	    llist = []
	    for g,seq in dict_seq.items():
		llist.append(seq[j])
		if seq[j] == "-":
		    break
	    if not "-" in llist:
		if len(list(set(llist))) == 2:
		    contig.list_SNP_pos.append(str(align_start+len(ref_seq[:j].replace("-", "")))) # convert to the physical pos in reference genome
		    for g,seq in dict_seq.items():
		        contig.dict_SNPs[g]+=seq[j]
		elif len(list(set(llist))) > 2:
		    contig.num_multi_locus_SNPs += 1
		contig.num_aligned_sites += 1

# list of genome names for writing phase file
list_genomes = []
with open("target_genomes.txt") as f, open("idfile.txt", "w") as f2:
    line = f.readline()
    while line:
        if not line.startswith("#") and len(line)>1:
	    list_genomes.append(line[:-1].split("\t")[0])
	    #list_genomes.append(line[:-1])
	    f2.write("%s" % line[1:].replace("\t", "_"))
        line = f.readline()


# write phase file for each contig
def writePhaseFile(ext):
  for name in list_contig_names:
    contig = dict_contigs[name]
    if len(contig.list_SNP_pos) > 1:
      with open("bialle_SNP.ref_%s.%s" % (name, ext), "w") as f:
        num_individuals = len(list_genomes)
	if "validation" in ext:
	    num_individuals += 1
        num_SNPs = len(contig.list_SNP_pos)
        str_SNP_pos = " ".join(contig.list_SNP_pos)
	f.write("%d\n%d\nP %s\n" % (num_individuals, num_SNPs, str_SNP_pos))
	for i in range(len(list_genomes)):
	    g = list_genomes[i]
	    f.write(contig.dict_SNPs[g] + "\n")
        if "validation" in ext:
	    f.write(contig.dict_SNPs[ref_genome]+"\n")
writePhaseFile("phase")

# write phase file for validation only, if ref not in targe_genomes
if ref_genome not in list_genomes:
    writePhaseFile("phase.validation")




# stat
summary = "contig_id\tcontig_size\tnum_sites_aligned(w_o_gap)\tnum/contig_size\tnum_SNPs\tnum/aligned\tnum_biallelic_SNPs\tnum/SNPs\tnum_alignment_blocks\tnum_individuals\n"
for name in list_contig_names:
    contig = dict_contigs[name]
    contig_length = contig.end - contig.start + 1 
    num_aligned_sites = contig.num_aligned_sites
    perc_n_aligned_sites = 100.0 * num_aligned_sites / contig_length

    num_SNPs = contig.num_multi_locus_SNPs + len(contig.list_SNP_pos)
    perc_n_snps = float("nan")
    if num_aligned_sites > 0:
        perc_n_snps = 100.0 * num_SNPs / num_aligned_sites

    num_biSNPs = len(contig.list_SNP_pos)
    perc_n_biSNPs = float("nan")
    if num_SNPs > 0:
	perc_n_biSNPs = 100.0 * num_biSNPs / num_SNPs
	
    num_align_blocks = len(contig.dict_alignments)
    num_individuals = len(contig.dict_SNPs)
	
    summary += "%s\t%d\t%d\t%.1f%%\t%d\t%.1f%%\t%d\t%.1f%%\t%d\t%d\n" % (name, contig_length, num_aligned_sites, perc_n_aligned_sites, \
		num_SNPs, perc_n_snps, num_biSNPs, perc_n_biSNPs, num_align_blocks, num_individuals)
summary = summary.replace("nan%", "--")
print summary

with open("log.STEP1.s3.SNP_summary.txt", "w") as f:
    f.write(summary)




