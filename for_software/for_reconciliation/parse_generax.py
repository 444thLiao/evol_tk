from os.path import *
from ete3 import Tree,TreeNode
from api_tools.itol_func import pie_chart
from bs4 import BeautifulSoup
import pandas as pd

def get_up_till_root(n):
    if n.up is not None:
        return [n.up] + get_up_till_root(n.up)
    else:
        return []
def get_trees_from_xml(soup):
    p = list(soup.find_all('phylogeny'))
    tree_dict = {}
    for _ in p:
        if _.parent.name =='spTree':
            tree_dict['speciesTree'] = _
        elif _.parent.name=='recGeneTree':
            tree_dict['reconciliatedGeneTree'] = _
    return tree_dict

def get_info(n):
    # if not 
    name = [_ for _ in n.children if _.name =='name'][0].text
    children = [_ for _ in n.children if _.name == 'clade']
    if len(children)==2:
        l,r = children
        return name,l,r
    else:
        return name,

def DFS_get_tree(root,par_node):
    results = get_info(root)
    par_node.name = results[0]
    if len(results)==1:
        # par node is a leaf, end
        
        return 
    elif len(results)==3:
        name,l,r = results
        l_node = TreeNode()
        r_node = TreeNode()
        par_node.add_child(l_node)
        par_node.add_child(r_node)
        return DFS_get_tree(l,l_node),DFS_get_tree(r,r_node)
    
def format_tree_into_ete(tree_dict):
    etetree_dict = {}
    for key,v in tree_dict.items():
        
        root = list(v.children)[1]
        tre = Tree()
        DFS_get_tree(root,tre) 
        etetree_dict[key] = tre
    return etetree_dict

def get_info_from_xml(xml_file):
    soup = BeautifulSoup(open(xml_file), 'xml')
    tree_dict = get_trees_from_xml(soup)
    etetree_dict = format_tree_into_ete(tree_dict)
    
    def format_dict(eventRec,name):
        # for eventsRec
        gene_node_xmlnode = eventRec.parent.previous_sibling.previous_sibling
        if name == 'speciesLocation':
            candidate_receptors = list(eventRec.parent.find_next_siblings('clade'))
            # within all sister Clade nodes, find the node with eventsRec and also the transferBack
            receptors = []
            for _ in candidate_receptors:
                child = list(_.findChildren('eventsRec',recursive=False))[0]
                # inhibit the recursive, otherwise, some sister Clade maybe contain a very 'depth' phylogeny which contain other eventsRec.
                if child.findChildren('transferBack'):
                    # only one of candidate_receptors is containing receptor
                    receptors.append(child.findChildren('transferBack')[0].attrs['destinationSpecies'])
                    
            assert len(receptors) == 1
            return gene_node_xmlnode.text,eventRec.attrs[name],receptors[0]
    # the donors and receptors should be the same.
    HGT_donors = {}
    for n in soup.find_all("branchingOut"):
        name_node,attr,receptor = format_dict(n,'speciesLocation')
        HGT_donors[(name_node,attr)] = (n,receptor)
    tranfer_pairs = [(k[1],v[1]) for k,v in HGT_donors.items()]
    return tree_dict,etetree_dict,HGT_donors,tranfer_pairs



def get_eventsize(xml_file):
    f = xml_file.replace("_reconciliated.xml","_speciesEventCounts.txt")
    df = pd.read_csv(f,sep=' ',header=None,index_col=0)
    df = df.drop(columns=5)
    df.columns = 'speciations,duplications,losses,transfers'.split(',')
    return df


def get_p2node(xml_file,stree=None,key='',return_all=False):
    p2node = {}
    p2node_transfer_receptor = {}
    if stree is None:
        tree_dict,etetree_dict,HGT_donors,tranfer_pairs = get_info_from_xml(xml_file)
        stree = etetree_dict['speciesTree']
    else:
        stree = Tree(stree,format=3)
    name2node = {} # children to parent
    for n in stree.traverse():
        if n.name:
            name2node[n.name] = n
                   
    soup = BeautifulSoup(open(xml_file), 'xml')
    all_speciations = soup.find_all("speciation") # speciation events
    all_s_nodes = [t.attrs['speciesLocation'] for t in all_speciations]
    remained_nodes = []
    for n in all_s_nodes:
        node = name2node[n]
        parents = [_.name for _ in get_up_till_root(node)]
        if set(parents).intersection(set(all_s_nodes)):
            continue
        else:
            remained_nodes.append(n)
    # get transfer events
    transfer_f = xml_file.replace('_reconciliated.xml','_transfers.txt')
    if getsize(transfer_f)!=0:
        transfer = [(row.split(' ')[0],row.split(' ')[1] )
                    for row in open(transfer_f).read().strip('\n').split('\n')]
        target_transfer = [_[1] for _ in transfer]
        speciation_but_gain_from_transfer = [_ for _ in remained_nodes if _ in target_transfer]
        remained_nodes = [_ for _ in remained_nodes if _ not in target_transfer]
        p2node[key] = remained_nodes
        p2node_transfer_receptor[key] = speciation_but_gain_from_transfer
    else:
        p2node[key] = remained_nodes
        p2node_transfer_receptor[key] = remained_nodes
    return p2node,p2node_transfer_receptor


def get_xml2itol(p2node,p2node_transfer_receptor,
                 ofile,
                 color_speciation='#ff0000',
                 color_transfer='#0000ff'):

    n2cat2v = {n:{"speciation":1} for n in p2node[p]}
    n2cat2v.update({n:{"transfer_receptor":1} for n in p2node_transfer_receptor[p]})
    text = pie_chart(n2cat2v,
                    cat2style={"speciation":color_speciation,
                                "transfer_receptor":color_transfer},
                    dataset_label="speciation")
    with open(ofile,'w') as f1:
        f1.write(text)

if __name__ == "__main__":
    gene = 'nxrA'
    odir = "./reconciliation/generax/focal_phyla"
    stree = "./reconciliation/generax/focal_phyla/stree.newick" 
    r_odir = join(odir,'nxrA_prot')
    xml_file = join(r_odir,f'reconciliations/{gene}_reconciliated.xml')
    
    
    sizedf = get_eventsize(xml_file)
    node2speciation_times = sizedf['speciations'].to_dict()
    from api_tools import pie_size_chart
    text = pie_size_chart(node2speciation_times,color='#e3772d')
    
    p2node,p2node_transfer_receptor = get_p2node(xml_file,stree,key='all')