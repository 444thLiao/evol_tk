from Bio import SeqIO



aln_file = "gene.prot.aln"
trimal_file = "gene.prot.trimal"

top1_r = next(SeqIO.parse(aln_file,format='fasta'))
top1_r_trimed = next(SeqIO.parse(trimal_file,format='fasta'))



