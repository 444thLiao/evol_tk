
import os,re
os.chdir('/mnt/storage1/jjtao/amplicon/F22FTSHMHT1079')

from subprocess import check_call
from Bio.Data import IUPACData
from Bio.Seq import Seq
from tqdm import tqdm
from collections import defaultdict
from Bio import SeqIO
from os.path import *
import copy
from ete3 import Tree
tmp_d = IUPACData.ambiguous_dna_values
tmp_d = {k.upper():v.upper() for k,v in tmp_d.items()}

def get_region(seq, f, r,include_primer=False):
    f_len, r_len = len(f), len(r)
    f = ''.join([l if len(tmp_d.get(l, l)) ==
                 1 else f"[{tmp_d.get(l,l)}]" for l in f])
    r = str(Seq(r).reverse_complement())
    f = ''.join([l if len(tmp_d.get(l, l)) ==
                 1 else f"[{tmp_d.get(l,l)}]" for l in f])
    r = ''.join([l if len(tmp_d.get(l, l)) ==
                 1 else f"[{tmp_d.get(l,l)}]" for l in r])
    if include_primer:
        matched_seq = re.search(f'({f}[ATCGatcg]*{r})', seq)
    else:
        matched_seq = re.search(f'{f}([ATCGatcg]*){r}', seq)
    if matched_seq is None:
        return None
    # return the sequence removing the primer
    return matched_seq.groups()[0]

def trim_aln_with_seq(ref_seq, regional_seq):
    # 0-coordinated, thus left closed and right closed
    aft2ori_idx = {}
    ori_i = aft_i = 0
    for s in ref_seq:
        if s not in 'actgnATCGN':
            ori_i += 1
        else:
            aft_i += 1
            ori_i += 1
        aft2ori_idx[aft_i] = ori_i
    ori_idx = list(range(0, len(ref_seq)))
    nogap_refseq = ''.join([_
                            for _ in ref_seq if _ in 'actgnATCGN'])
    aft_idx = list(range(0, len(nogap_refseq)))

    idx = nogap_refseq.find(regional_seq)
    if idx == -1:
        raise IOError
    start = aft2ori_idx[idx]
    end = aft2ori_idx[idx+len(regional_seq)]
    return start, end

def get_matched(pre_sequences, primer_set,ofile):
    matched_ = defaultdict(int)
    for r in pre_sequences:
        pos = [_ for _,bp in enumerate(r.seq) if bp !='-']
        new_s = [bp for _,bp in enumerate(r.seq) if bp !='-']
        pos_d = dict(zip(range(len(new_s)),pos))
        region = get_region(''.join(new_s).upper(),
                            primer_set['f'],
                            primer_set['r'],
                            include_primer=False)
        
        if not region:
        # print(r.id)
            continue
        s, e = trim_aln_with_seq(''.join(new_s).upper(), region)
        ns = pos_d[s]
        ne = pos_d[e]
        matched_[(ns,ne)] +=1
    ### manual check whether it has only one start & end. Otherwise it might have some problems
    if len(matched_)!=1:
        print(matched_,primer_set)
        return
    s,e = list(matched_.keys())[0]    
    trimmed_seqs = []
    for r in pre_sequences:
        sub_r = copy.deepcopy(r)[s:e]
        #sub_r.seq = Seq(''.join([_ for _ in sub_r.seq if _!='-']))
        if len([_ for _ in sub_r.seq if _!='-'])<=50:
            continue
        trimmed_seqs.append(sub_r)    
    with open(ofile,'w') as f:
        SeqIO.write(trimmed_seqs,f,'fasta-2line')

def batch_iter(iter, batch_size):
    # generating batch according batch_size
    iter = list(iter)
    n_iter = []
    batch_d = 0
    for batch_u in range(0, len(iter), batch_size):
        if batch_u != 0:
            n_iter.append(iter[batch_d:batch_u])
        batch_d = batch_u
    n_iter.append(iter[batch_d : len(iter) + 1])
    return n_iter




primer_set = {"f": "TGCGAYCCSAARGCBGACTC", #second
              "r": "ATSGCCATCATYTCRCCGGA"}
'''
PolF 	TGCGAYCCSAARGCBGACTC
PolR	ATSGCCATCATYTCRCCGGA
BH333f	GARGANGTBATGAAGGTCGG
BH773r	TGCTGCACGATRTTGTCGC
BRB2106f	CCGRTSACGCCBGACAAG
BRB2516r	TGTCGCCCTTCYTGACGAYR
338F 	ACTCCTACGGGAGGCAGCA
806R	GGACTACHVGGGTWTCTAAT
'''

nifH_full_aln = '/home-user/jjtao/Rhizobiales2/gene_tree/nt-nifH/full_length2/nifH.aline'
pre_sequences = list(SeqIO.parse(nifH_full_aln,'fasta'))

cmds = []
fname2primerset = {"F22FTSHMHT1079_337f773r":{"f": "GARGANGTBATGAAGGTCGG", 
                                              "r": "TGCTGCACGATRTTGTCGC",
                                              'seqs':pre_sequences},
                   "F22FTSHMHT107901_polFR":{"f": "TGCGAYCCSAARGCBGACTC",
                                             "r": "ATSGCCATCATYTCRCCGGA",
                                             'seqs':pre_sequences},}

for fname,primer_set in fname2primerset.items():
    pre_sequences = primer_set['seqs']
    ofile_name = f'/mnt/storage1/jjtao/amplicon/F22FTSHMHT1079/{fname}/output/dada2_output'
    get_matched(pre_sequences, {k:v for k,v in primer_set.items() if k in 'fr'}
                ,ofile_name+'/trimmed_ref_nirH.aln')

    ref = list(SeqIO.parse(ofile_name+'/trimmed_ref_nirH.aln','fasta'))
    rep_fa = ofile_name+'/rep.fa'
    all_asv = list(SeqIO.parse(rep_fa,'fasta'))
    
    ref = [r for r in ref if len([_ for _ in r.seq if _!='-'])>50]
    ref_phy = ofile_name+'/tmp.phy'
    with open(ref_phy,'w') as f1:
        f1.write(f"{len(ref)} {len(ref[0].seq)}\n")
        for _r in ref:
            f1.write(f"{_r.id}{' '*10}{_r.seq}\n")  
    with open(ref_phy.replace('.phy','.aln'),'w') as f1:
        for _r in ref:
            f1.write(f">{_r.id}\n{_r.seq}\n")

    cmd2 = f"/home-user/thliao/anaconda3/bin/FastTree {ref_phy.replace('.phy','.aln')} > {ref_phy.replace('.phy','.newick')}"
    os.system(cmd2)            
    
    gene_tre = Tree(ref_phy.replace('.phy','.newick'))
    gene_tre.resolve_polytomy()
    c = 1
    for n in gene_tre.traverse():
        if not n.is_leaf():
            n.name = f"INode{c}"
            c+=1
    gene_tre.write(outfile=ref_phy.replace('.phy','_renamed.newick'),format=3)
    if exists(join(ofile_name,'papara_log.aln')):
        os.system(f"rm -r {join(ofile_name,'papara_log.aln')}")     
    if not exists(f"{ofile_name}/query.fasta"):
        cmd1 = f"cd {ofile_name} && /mnt/home-user/thliao/software/papara_nt/papara -t {ref_phy.replace('.phy','_renamed.newick')} -s {ref_phy} -q {rep_fa} -r -n aln -j 30 "
        check_call(cmd1,shell=1)
        cmd2 = f"cd {ofile_name} && /mnt/home-user/thliao/anaconda3/bin/epa-ng --split {ref_phy} papara_alignment.aln --redo"
        check_call(cmd2,shell=1)
    
    used_asv_ids = [_.id for _ in all_asv]
    ref_ids = [_.id for _ in ref]
    all_aln = {_.id:_ for _ in SeqIO.parse(f"{ofile_name}/query.fasta",'fasta')}
    all_aln.update({_.id:_ for _ in SeqIO.parse(f"{ofile_name}/reference.fasta",'fasta')})
    n_iter = batch_iter(used_asv_ids,len(used_asv_ids)//50)
    for idx,n in tqdm(enumerate(n_iter),total=len(n_iter)):
        aln_subset = [all_aln[asv] for asv in n] + [all_aln[_] for _ in ref_ids]
        if not exists(f"{ofile_name}/asv_tree"):
            os.makedirs(f"{ofile_name}/asv_tree")        
        with open(f'{ofile_name}/asv_tree/ASV_{idx}.aln','w') as f:
            SeqIO.write(aln_subset,f,'fasta-2line')
        cmd2 = f"/home-user/thliao/anaconda3/bin/FastTree {ofile_name}/asv_tree/ASV_{idx}.aln > {ofile_name}/asv_tree/ASV_{idx}.tree"
        cmds.append(cmd2)
        #os.system(cmd2)


