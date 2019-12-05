from self_substitution_calculator import get_codon,translate

base4 = {'A','C','T','G'}
def generate_point_muta_all_possible(seq):
    assert ' ' not in seq
    seq = seq.upper()
    codons = get_codon(seq)
    
    others = []
    for idx,base in enumerate(seq):
        other_possible = base4.difference(set(base))
        for ob in other_possible:
            new_seq = seq[:idx] + ob +seq[idx+1:]
            others.append(new_seq)
    return others
        
def calculate_num_sites(seq):
    assert ' ' not in seq
    # original_protein_seq = translate(seq)
    codons = get_codon(seq)
    total_syn = 0
    total_nonsyn = 0
    for c in codons:
        ori = translate(c)
        all_possible_seq = generate_point_muta_all_possible(c)
        protein_seqs = [translate(_)
                   for _ in all_possible_seq]
        syn = len([_ for _ in protein_seqs if _ == ori])
        nonsyn = len(protein_seqs) - syn
        total_syn += syn/3
        total_nonsyn += nonsyn/3
    return total_syn,total_nonsyn

def get_pair_num_sites(seq1,seq2):
    total_syn1,total_nonsyn1 = calculate_num_sites(seq1)
    total_syn2,total_nonsyn2 = calculate_num_sites(seq2)
    
    final_syn = (total_syn1+total_syn2)/2
    final_nonsyn = (total_nonsyn1+total_nonsyn2)/2
    return final_syn,final_nonsyn

def get_pair_num_sites_with_freq(seq1,seq2):
    total_syn1,total_nonsyn1 = calculate_num_sites(seq1)
    total_syn2,total_nonsyn2 = calculate_num_sites(seq2)
    
    final_syn = (total_syn1+total_syn2)/2
    final_nonsyn = (total_nonsyn1+total_nonsyn2)/2
    return final_syn,final_nonsyn

def get_3x4(seq):
    # useless for now
    codons = get_codon(seq)
    pos_total = {pos+1:{base:0 for base in base4} 
                for pos in range(0,3)
                }
    for codon in codons:
        for pos,base in enumerate(codon):
            pos_total[pos+1][base] +=1
    pos_freq = {pos:{base: val/sum(pos_total[pos].values()) for base,val in pos_total[pos].items()} 
                for pos in pos_total}
    return pos_freq
    
seq1 = 'TCAACTGAGATGTGTTTA'
final_syn,final_nonsyn = get_pair_num_sites(seq1,seq1)
final_syn,final_nonsyn
seq1 = 'TCAACTGAGATGTGTTTA'
seq2 = 'TCAACAGAGATATGTCTA'
# different to yn00
final_syn,final_nonsyn = get_pair_num_sites(seq1,seq2)
final_syn,final_nonsyn

seq1 = 'AGTACTGAGATGTGTTTA'
seq2 = 'TCTACAGAGATATGTCTA'
# different to yn00
final_syn,final_nonsyn = get_pair_num_sites(seq1,seq2)
final_syn,final_nonsyn

seq1 = 'AGTACTGAGATGTGTTTA'
seq2 = 'AGTACTGAGATGTGTTTA'
# different to yn00
final_syn,final_nonsyn = get_pair_num_sites(seq1,seq2)
final_syn,final_nonsyn


seq1 = 'ATGAAACCCGGGTTT'
seq2 = 'ATGAAACCCGGGTTT'
# same as yn00
final_syn,final_nonsyn = get_pair_num_sites(seq1,seq2)
final_syn,final_nonsyn
