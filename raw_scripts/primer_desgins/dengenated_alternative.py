from Bio.Data import IUPACData
dgbase = IUPACData.ambiguous_dna_values

def expand_dg_primer(inseq):
    record_ = [[]]
    for _,b in enumerate(inseq):
        if b in 'actg':
            for base_seq in record_:
                base_seq.append(b)
        else:
            seed_seqs = record_[::]
            record_ = []
            for base_seq in seed_seqs:
                for dg_b in dgbase[b]:
                    record_.append(base_seq+[dg_b])
    ep_primers = [''.join(_) for _ in record_]
    return ep_primers

target_seq = "CGATCAGYTTGGAYTTCTGRACCTG"
a = expand_dg_primer(target_seq)

