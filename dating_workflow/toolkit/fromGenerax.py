from os.path import *
from ete3 import Tree
from for_software.for_reconciliation.parse_generax import get_p2node
from dating_workflow.step_script.quick_sampling import *
from api_tools import to_binary_shape,to_color_range,pie_chart

def generax2mcmctree(xml_file,stree,gene,dating_o,
                     calbration_file,
                     genome2cog25={}
                     ):
    """
    generate files for mcmctree, including 
    1. list of used genomes
    2. used genomes (itol annotation)
    3. constructed species tree topology for dating
    4. target genomes in the complete phylogeny of phylum (itol annotation)
    5. calibration file and tree file with calibrations information
    """
    # xml_file = join(r_odir,f'reconciliations/{gene}_reconciliated.xml')
    # stree = f"./trees/iqtree/{phylum_name}.reroot.newick" 
    
    phylum_name = xml_file.split('/')[-4]
    st =  Tree(stree,format=3)
    tmp_name = phylum_name+'_'+gene
    _p2node,_p2node_transfer_receptor = get_p2node(xml_file,stree,key=tmp_name)
    target_nodes = list(_p2node.values())[0] + list(_p2node_transfer_receptor.values())[0]
    
    must_in_genomes = open("/mnt/home-backup/thliao/cyano_basal/rawdata/assembly_ids.list").read().strip('\n').split('\n')
    # new calibrations are /mnt/home-backup/thliao/cyano/ref_genomes_list.txt
    cluster2genomes = get_cluster(stree.replace('.reroot.newick','.clusterd.list'))
    g2cluster = {v:c for c,d in cluster2genomes.items() for v in d}
    retained_ids = sampling(st,target_nodes,must_in = must_in_genomes,node2cluster=g2cluster,genome2cog25=genome2cog25)

    text = to_binary_shape({g:['keep'] for g in retained_ids},
                        {"keep":{"color":"#88b719"}})
    text = to_color_range({g:'keep' for g in retained_ids},
                          {"keep":"#88b719"} )
    with open(join(dating_o,f'id_list/{phylum_name}_{gene}.txt'),'w') as f1:
        f1.write(text)
    with open(join(dating_o,f'id_list/{phylum_name}_{gene}.list'),'w') as f1:
        f1.write('\n'.join(retained_ids))
    print(phylum_name,
          len(st.get_leaf_names()),
          len(retained_ids))
    st.copy()
    st.prune(retained_ids)
    with open(join(dating_o,f'species_trees/{phylum_name}_{gene}.newick'),'w') as f1:
        f1.write(st.write(format=9))
    
    # draw target nodes
    LCA_nodes =[]
    for name in target_nodes:
        n = [n for n in st.traverse() if n.name == name][0]
        l1 = n.children[0].get_leaf_names()[0]
        l2 = n.children[1].get_leaf_names()[0]
        LCA_nodes.append(f"{l1}|{l2}")
        
    text = pie_chart({n:{'speciation':1} for n in LCA_nodes},
                           {"speciation":"#ff0000"},
                           dataset_label='GeneRax results')
    with open(join(dating_o,f'target_nodes/{phylum_name}_{gene}.txt'),'w') as f1:
        f1.write(text)

    # new set file set14
    # set14_f = './dating/calibration_sets/scheme1/cal_set14.txt'
    
    c = 'GCA_000011385.1' 
    n = [_ for _ in st.children if c  not in _.get_leaf_names()][0]
    final_text = open(calbration_file).read().replace('GCA_002239005.1',n.get_leaf_names()[0])
    with open(join(dating_o,f'calibrations/{phylum_name}_{gene}_set14.txt'),'w') as f1:
        f1.write(final_text)
        