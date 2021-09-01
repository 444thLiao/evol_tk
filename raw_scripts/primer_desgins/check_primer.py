from collections import defaultdict
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq
import re

## alternative choice
f_primer = 'CGCGAGAACAAGGCGCAGGA'
r_primer = "TCCAGCTCCTTGTCGGTCTG"

## best choice
f_primer = "CAGACCGACAAGGAGCTGGA"
r_primer = "TTCATGATGCCGTGCTCCATCA"


f_primer = "CGCGTTTTACGGCAAGGGC"
r_primer = "GCCTTGGCGTGCAGGATCAG"
r_primer_rc = str(Seq(r_primer).reverse_complement())


records = list(SeqIO.parse('ref_up.aln','fasta'))
anno = open('./seq_annotations.txt').read().split('\n')[1:-1]
seq2anno = {row.split('\t')[0]:row.split('\t')[1] for row in anno}
type2seq = defaultdict(list)
for _seq,a in seq2anno.items():
    type2seq[a].append(_seq)

def get_match_count(f_primer,r_primer):
    r_primer_rc = str(Seq(r_primer).reverse_complement())
    match_count = defaultdict(int)
    for r in records:   
        seq = ''.join([_ for _ in str(r.seq) if _!='-'])
        region = re.findall(f"{f_primer}.+{r_primer_rc}",seq)
        if region:
            match_count[seq2anno[r.id]]+=1
    return match_count

target_seqs = []
match_count = defaultdict(int)
for r in records:   
    seq = ''.join([_ for _ in str(r.seq) if _!='-'])
    region = re.findall(f"{f_primer}.+{r_primer_rc}",seq)
    if region:
        match_count[seq2anno[r.id]]+=1
    target_seqs.append(r.id)

from skbio import DNA, TabularMSA
msa = TabularMSA([DNA(str(_.seq)) for _ in records])
msa.index = [_.id for _ in records]

sym_conserve = msa.loc[type2seq['Sym'],:].conservation()
fl_conserve = msa.loc[type2seq['FL'],:].conservation()


valid_pos = []
for pos,v in enumerate(fl_conserve):
    if v == 1:
        aa = str(msa.loc[type2seq['FL'][0],pos])
        _v = sym_conserve[pos]
        if _v <1 and str(_v)!='nan':
            freq = str(msa.loc[type2seq['Sym'],pos]).count(aa)/len(type2seq['Sym'])
            valid_pos.append((pos,aa,freq))
            
            
# scan search
from tqdm import tqdm
length = 20

alternative_primer = []
for start_pos in tqdm(range(msa.shape[1]-19)):
    end_pos = start_pos + length
    sub_msa = msa.iloc[:,start_pos:end_pos]
    
    sym_cons_seq = sub_msa.loc[type2seq['Sym'],:].consensus()
    fl_cons_seq = sub_msa.loc[type2seq['FL'],:].consensus()
    if str(sym_cons_seq)==str(fl_cons_seq):
        # consensus sequences within each group are different
        continue
    if fl_cons_seq.has_gaps():
        # no gap
        continue
    if fl_cons_seq.gc_content() <=0.4 or fl_cons_seq.gc_content() >=0.6:
        continue
    diff_sites = [_1 for _1,_2 in zip(str(fl_cons_seq),str(sym_cons_seq))
                  if _1!=_2]
    num_diff_sites = len(diff_sites)
    alternative_primer.append((start_pos,
                               str(fl_cons_seq),
                               str(sym_cons_seq),
                               num_diff_sites))
    
alternative_primer = sorted(alternative_primer,key=lambda x:x[-1],reverse=True)

get_match_count('GCACAAGATCCACAACAACG',"CGGAATTGGCGTATTTCAGA")





# alternative primer

from skbio import DNA, TabularMSA
from Bio import SeqIO
records = list(SeqIO.parse('./related_rpoB.aline','fasta'))
msa = TabularMSA([DNA(str(_.seq).upper()) for _ in records])
msa.index = [_.id for _ in records]

g1 = [_.id for _ in msa if 'Afipia' in _.id]
g2 = [_.id for _ in msa if 'Bradyrhizobium' in _.id]
g3 = [_.id for _ in msa if 'Rhodopseudomonas' in _.id]

p2c = {}
length = 20
for start_pos in tqdm(range(msa.shape[1]-19)):
    end_pos = start_pos + length
    sub_msa = msa.iloc[:,start_pos:end_pos]
    
    g1_cons_seq = sub_msa.loc[g1,:].consensus()
    g2_cons_seq = sub_msa.loc[g2,:].consensus()    
    g3_cons_seq = sub_msa.loc[g3,:].consensus() 
    if sub_msa.loc[g2,:].conservation().min() >=0.85:
        p2c[start_pos] = [str(g1_cons_seq),str(g2_cons_seq),str(g3_cons_seq)]
    if g2_cons_seq.has_degenerates():
        print(g2_cons_seq)

