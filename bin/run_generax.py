#!/usr/bin/env python
"""
run generax using reference tree and a gene tree
"""
import sys
from os.path import dirname
sys.path.insert(0, dirname(__file__))
import click
from api_tools.for_tree.format_tree import renamed_tree, root_tree_with, add_cal_api, Tree, read_tree, sort_tree, earse_name, draw_cal_itol
import os
from api_tools.itol_func import to_node_symbol
from os.path import exists,join,basename
from dating_workflow.step_script import process_path

infaa = ''



collect_parameters = []
redo = True
for gene in genes:
    record_list = list(SeqIO.parse(join(base_odir,'..',f"{gene}.faa"),'fasta'))
    l2g_list = join(generax_odir,f'{gene}_locus2genome.txt')
    if not exists(dirname(l2g_list)):
        os.system(f'mkdir -p {dirname(l2g_list)}')
    with open(l2g_list,'w') as f1:
        f1.write('\n'.join([f"{r.id}\t{r.id.split('_')[0]}" for r in record_list]))
    
    ofaa = join(generax_odir,f"{gene}.faa")
    for r in record_list:
        r.name=r.description = ''
    sub_record_list = [_ for _ in record_list if _.id.split('_')[0] in all_genomes]
    SeqIO.write(sub_record_list,
                open(ofaa,'w'),
                format='fasta-2line')
    aln_f = ofaa.replace('.faa','.aln')
    if not exists(aln_f):
        os.system(f"ginsi {ofaa} > {aln_f}")
    gtree_f = join(base_odir.replace('/generax','/') ,gene+'.newick')
    new_gtree = join(base_odir,gene+'.newick')
    new_gt = Tree(gtree_f)
    new_gt.resolve_polytomy()
    with open(new_gtree,'w') as f1:
        f1.write(new_gt.write())
    collect_parameters.append((gene,new_gtree,aln_f,l2g_list))
        
with open(join(generax_odir,f"genrax.ctl"),'w') as f1:
    f1.write(f"""[FAMILIES] # this is a comment""" + '\n')
    for gene,sgtree,align,mapping in collect_parameters:
        f1.write(f"""
- {gene}
starting_gene_tree = {sgtree}
alignment = {align}
mapping = {mapping}
subst_model = LG+I+G4    """)
        
# from bin.format_newick import renamed_tree
# t = Tree(stree,3)
# t.set_outgroup(t.get_common_ancestor(["GCA_011774645.1","GCA_014533965.1"]))
#t.set_outgroup(t.get_common_ancestor(['GCA_003242895.1','GCA_000011385.1']))  # topo1
with open(join(generax_odir,'stree.newick'),'w') as f1:
    f1.write(genome_tree.write(format=3))
cmd = f"mpiexec -np 32 generax -s {join(generax_odir,'stree.newick')} -p {join(generax_odir,'output')} -f {join(generax_odir,f'genrax.ctl')} --seed 100 --unrooted-gene-tree --per-family-rates"
print(cmd)