
import os
os.chdir('/mnt/storage1/jjtao/amplicon/F22FTSHMHT1079')

from tqdm import tqdm
from collections import defaultdict
from Bio import SeqIO
from os.path import *
from ete3 import Tree
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
from collections import Counter
from glob import glob


def find_LCAs(leaves,tre,ref_ids):
    tre = tre.copy()
    lca = tre.get_common_ancestor(leaves)
    new_ids = [_ for _ in lca.get_leaf_names() if _ in ref_ids]
    larger_one = sorted([set(new_ids),set(leaves)],key=lambda x:len(x))
    if len(set(new_ids)&set(leaves))/len(larger_one[-1])<=0.95:
        lca_list = []
        iter_nodes = list(lca.children)
        while iter_nodes:
            for n in list(iter_nodes):
                new_ids = [_ for _ in n.get_leaf_names() if _  in ref_ids]
                if len(new_ids) == 0:
                    iter_nodes.remove(n)
                    continue
                larger_one = sorted([set(new_ids),set(leaves)],key=lambda x:len(x))
                if len(set(leaves)&set(new_ids))/len(larger_one[-1])<=0.95:
                    iter_nodes.remove(n)
                    iter_nodes += list(n.children)
                else:
                    iter_nodes.remove(n)
                    lca_list.append(n)
        return lca_list
    else:
        return [lca]

fname2primerset = {
                   "F22FTSHMHT1079_337f773r":{"f": "GARGANGTBATGAAGGTCGG",
                                              "r": "TGCTGCACGATRTTGTCGC"},
                #    "F22FTSHMHT107901_polFR":{"f": "TGCGAYCCSAARGCBGACTC",
                #                              "r": "ATSGCCATCATYTCRCCGGA"},
                   }

for fname,primer_set in fname2primerset.items():
    ofile_name = f'/mnt/storage1/jjtao/amplicon/F22FTSHMHT1079/{fname}/output/dada2_output'
    reference_names = list(SeqIO.parse(ofile_name+'/trimmed_ref_nirH.aln','fasta'))
    reference_names = [r.id for r in reference_names if len([_ for _ in r.seq if _!='-'])>50]

    name2node = {_.name:_ for _ in Tree(ofile_name+'/tmp_renamed.newick',format=3).traverse()}
    group2IN = defaultdict(list)
    for row in open(f'{ofile_name}/node2color').read().strip().split('\n'):
        group2IN[row.split('\t')[1]].append(name2node[row.split('\t')[0]])

    asv2groups = {}
    for tpath in tqdm(glob(f"{ofile_name}/asv_tree/ASV_*.tree")):
        tre = Tree(tpath)
        asv_l = [_.name for _ in tre.get_leaves() if _.name not in reference_names]
        for g,nodes in group2IN.items():
            for n in nodes:
                ref_gids = n.get_leaf_names()
                if len(ref_gids)==1:
                    continue
                lca_list = find_LCAs(ref_gids,tre,reference_names)
                for lca in lca_list:
                    for l in lca.get_leaf_names():
                        if l not in reference_names:
                            asv2groups[l] = g
        remaining_asv = set(asv_l).difference(set(asv2groups.keys()))
        for a in remaining_asv:
            asv2groups[a] = 'Unassigned'
    with open(f'{ofile_name}/assigned_asv.txt','w') as f1:
        for asv,group in asv2groups.items():
            f1.write(f'{asv}\t{group}\n')
    print(Counter(asv2groups.values()))
