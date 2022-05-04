from collections import defaultdict
from ete3 import Tree
from glob import glob
from os.path import *

##! read the branch estimated by codeml/baseml deposited in the 'in.BV'
def read_phy(phy_file,names=[]):
    rows = [_ for _ in open(phy_file).read().split('\n') if _]
    _d = defaultdict(dict)
    idx = -1
    for row in rows:
        a,b = [_ for _ in row.split(' ') if _]
        try:
            NoGenome = int(a)
            Nosites = int(b)
            idx+=1
        except:
            genome = a
            seq = b
            if set(seq) == set('-'):
                continue
            _d[idx][genome] = seq
    if names:
        _d = {names[idx]:v for idx,v in _d.items()}
    return _d

def get_branches_dist(t,leaves=[]):
    tre = Tree(t)
    if leaves:
        tre.prune(leaves)
    dist = tre.get_farthest_leaf()[1]
    return dist,tre

def read_branches_dist(rst2,leaves=[]):
    contexts = [_ for _ in open(rst2).read().split('\n') if _]
    all_tre_texts = [idx for idx,_ in enumerate(contexts) if _ == contexts[0]]
    all_tre_texts = [contexts[idx+1] for idx in all_tre_texts]
    assert len(all_tre_texts) == 1
    dist,tre = get_branches_dist(all_tre_texts[0],leaves=leaves)
    return dist,tre

def read_partition(partition_file):
    rows = [_ for _ in open(partition_file).read().split('\n') if _]
    name2length = {}
    names = []
    for r in rows:
        vals = [_ for _ in r.split(' ') if _ not in [',','=','']]
        gene = vals[1]
        s,e = vals[-1].split('-'); s,e = int(s),int(e)
        length = (e-s)+1
        name2length[gene] = length
        names.append(gene)
    return names,name2length

def get_rates(indir,write_tree=''):
    
    target_dir = dirname(indir)+'/tmp_files'
    ctl = target_dir+'/01_mcmctree_modify.ctl'
    row = [_ for _ in open(ctl).read().split('\n') if _.startswith('seqfile = ')][0]
    phy_file = row.strip().split(' ')[-1]
    partition_file = phy_file.replace('.phy','.partition')
    genes,name2length = read_partition(partition_file)
    gene2genome2seq = read_phy(phy_file,genes)
    _gene2rates = {}
    _gene2s = {}
    for _ in glob(f"{target_dir}/tmp*/tmp*.log"):
        name = _.split('/')[-2]
        row = [_ for _ in open(_).read().strip().split('\n') if 'ns = ' in _ and 'ls = ' in _][0]
        ns,sites,patterns = int(row.split('\t')[0].strip().split(' ')[-1]),int(row.split('\t')[-1].split(' ')[-1]),0
        gene = genes[int(name.replace('tmp',''))-1]
        _gene2s[gene] = (ns,sites,patterns)
        dist,tre = read_branches_dist(dirname(_)+'/rst2',list(gene2genome2seq[gene]))
        _gene2rates[gene] = dist
        if write_tree:
            tre.write(outfile=write_tree+f"/{gene}.newick")
    return _gene2s,_gene2rates

# ! 2. 



if __name__ == '__main__':
    # example of usage. Major functions
    # ! 1. read the branch length (represent the number of substituions per site. )
    name2s,name2rates = get_rates('',write_tree='outtree')
    
    # ! 2.
    
    