from ete3 import Tree,NCBITaxa
from Bio import Entrez
import io
from tqdm import tqdm
from os.path import *
from subprocess import check_call
ncbi = NCBITaxa()
def reformat(s):
    a = s.split('_')[-1]
    if not '_' in s:
        return s
    try:
        float(a)
        return s
    except:
        if len(s.rpartition('_')[-1]) == 1:
            return s
        else:
            return s.rpartition('_')[0]
    




for tree_file in ['nr_retrieve_hao/filtered_by_kegg.faa_aln.dir/iqtree.no_trim.treefile/K10535.sorted.newick',
                  'nr_retrieve_nxrA/cluster_95_filtered_lengths.fa_aln.dir/iqtree.treefile/K00370.sorted.newick',
                  'nr_retrieve_amoB/filtered_by_kegg.faa_aln.dir/iqtree.treefile/K10945.sorted.newick',
                  'nr_retrieve_amoC/filtered_by_kegg.faa_aln.dir/iqtree.treefile/K10946.sorted.newick']:
    basedir = '/home-user/thliao/project/nitrogen_cycle/fetch_genes/query_result'
    tree_file = join(basedir,tree_file)
    if not exists(tree_file):
        print(tree_file)
        
        #continue

    t = Tree(tree_file,format=1)
    all_ids = t.get_leaf_names()
    all_ids = [reformat(_) for _ in all_ids]
    basedir = dirname(tree_file)
    with open(join(basedir,'used_ids.list'),'w') as f1:
        f1.write('\n'.join(all_ids))

    infile = join(basedir,'used_ids.list')
    check_call(f'python3 /home-user/thliao/script/evolution_relative/global_search/get_info.py -i {infile} -o {dirname(infile)}',shell=1)

