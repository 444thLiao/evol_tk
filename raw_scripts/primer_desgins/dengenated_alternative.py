from itertools import product
dgbase = {
'A': ['A'],
'C': ['C'],
'G': ['G'],
'T': ['T'],
'R': ['A', 'G'],
'Y': ['C', 'T'],
'S': ['G', 'C'],
'W': ['A', 'T'],
'K': ['G', 'T'],
'M': ['A', 'C'],
'B': ['C', 'G', 'T'],
'D': ['A', 'G', 'T'],
'H': ['A', 'C', 'T'],
'V': ['A', 'C', 'G'],
'N': ['A', 'C', 'G', 'T'],
'I': ['A', 'C', 'G', 'T'],}

def expand_dg_primer(seq):
   """return a list of all possible k-mers given a degenerate base"""
   return list(map("".join, product(*map(dgbase.get, seq))))

target_seq = "CGATCAGYTTGGAYTTCTGRACCTG"
a = expand_dg_primer(target_seq)

