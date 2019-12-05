aa_dict = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S",
           "TCA": "S", "TCG": "S", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
           "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CTT": "L", "CTC": "L",
           "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
           "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R",
           "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
           "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N",
           "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
           "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A",
           "GCA": "A", "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E",
           "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}


def translate(seq,offset=0):
    """Does the actual translating"""
    tmp_str = ''
    for x in range(int(offset), len(seq), 3):
        codon = seq[x:x + 3]
        if len(codon) < 3:
            break
        elif "N" in codon:
            tmp_str += "X"
        else:
            tmp_str += "%s" % aa_dict.get(codon,'X')
    return tmp_str


def get_codon(seq,offset=0):
    """Does the actual translating"""
    codons = []
    for x in range(int(offset), len(seq), 3):
        codon = seq[x:x + 3]
        
        if len(codon) < 3:
            break
        codons.append(codon)
    return codons