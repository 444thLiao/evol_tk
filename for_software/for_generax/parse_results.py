from os.path import *
from ete3 import Tree
from api_tools.itol_func import pie_chart
from bs4 import BeautifulSoup

def get_up_till_root(n):
    if n.up is not None:
        return [n.up] + get_up_till_root(n.up)
    else:
        return []

def get_tree_from_xml(xml_file):
    soup = BeautifulSoup(open(xml_file), 'xml')
    all_speciations = soup.find_all("speciation") # speciation events
    all_s_nodes = [t.attrs['speciesLocation'] for t in all_speciations]
    # construct tree from xml_file (more exact)
    def get_info(n):
        # if not 
        name = [_ for _ in n.children if _.name =='name'][0].text
        children = [_ for _ in n.children if _.name == 'clade']
        if len(children)==2:
            l,r = children
            return name,l,r
        else:
            return name,
    from ete3 import TreeNode 
    p = list(soup.find_all('phylogeny'))
    p = [_ for _ in p if not _.attrs][0]
    root = list(p.children)[1]
    tre = Tree()
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
        # else:
        #     import pdb;pdb.set_trace()
        #     print(len(results))
        #     return
    DFS_get_tree(root,tre)      
    return tre

def get_p2node(xml_file,stree,key=''):
    p2node = {}
    p2node_transfer_receptor = {}
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
    p2node,p2node_transfer_receptor = get_p2node(xml_file,stree,key='all')