import pandas as pd
from os.path import *
import os
from Bio import SeqIO
import io
from subprocess import check_call
from ete3 import NCBITaxa
import plotly.express as px
from Bio import Entrez
from api_tools.itol_func import * 
from global_search.thirty_party.EntrezDownloader import EntrezDownloader

ncbi = NCBITaxa()

odir = './nr_retrieve_hao'
infile = './nr_retrieve_hao/protein2INFO.xlsx'
kofam_scan = '/home-user/thliao/software/kofamscan/exec_annotation'
ko = 'K10535'
protein_df = pd.read_excel(infile,index_col=0)
tree_exe = 'iqtree'


odir = './nr_retrieve_amoC'
infile = './nr_retrieve_amoC/protein2INFO.tab'
kofam_scan = '/home-user/thliao/software/kofamscan/exec_annotation'
ko = 'K10946'
protein_df = pd.read_csv(infile,sep='\t',index_col=0)
tree_exe = 'iqtree'
filter_id_txt = join(odir,'removed_ids.txt')


def run(cmd):
    check_call(cmd,shell=True)
    
def add_ref_seq(sub_df,used_records):
    used_ids = [str(_.id) for _ in used_records]
    id2info = {}
    record_need_dropped_ids = []
    for _,row in sub_df.iterrows():
        aa_id = row['AA accession'].strip()
        gene_name = row['gene name']
        seq = row['seq']
        info = row['phylum/class']
        if aa_id in used_ids:
            record_need_dropped_ids.append(aa_id)
            
        used_records.append(SeqIO.read(io.StringIO(f'>{aa_id}_{gene_name}\n{seq}'),format='fasta'))
        id2info[f'{aa_id}_{gene_name}'] = info
    final_records = [_ for _ in used_records if str(_.id) not in record_need_dropped_ids]    
    return final_records,id2info

color_scheme = {'type':{'NOB': '#e41a1c', 'comammox': '#edc31d', 
                        'AOB': '#bad5b9', 'AOA': '#358f0f'},
                'phylum/class':{'Thaumarchaeota': '#358f0f',
                                'Nitrospirae': '#edc31d',
                                'Gammaproteobacteria': '#78fce0',
                                'Chloroflexi': '#e41a1c',
                                'Betaproteobacteria': '#956cb4',
                                'Alphaproteobacteria': '#8c613c',
                                'Actinobacteria':'#11FF11',
                                'Planctomycetes':'#FF66bb',
                                }}

def get_colors_general(ID2infos,now_info2style={}):
    _c = color_scheme.copy()
    _c = _c['phylum/class']
    colors = px.colors.qualitative.Dark24 + px.colors.qualitative.Light24
    remained_colors = [c for c in colors if c not in now_info2style.values()]
    info2style = {}
    for v in set(ID2infos.values()):
        if v in _c:
            info2style.update({v:_c[v]})
        else:
            one_color = remained_colors.pop(0)
            info2style.update({v:one_color})
    return ID2infos,info2style

def write2colorbranch_clade(id2info,odir,info2color,treefile, unique_id,info_name='type',
                            **kwargs):
    content = to_color_Clade(id2info,info2color,treefile,info_name,**kwargs)
    info_name = info_name.replace('/','_')
    with open(join(odir, f'{unique_id}_{info_name}_colorbranch_clade.txt'), 'w') as f1:
        f1.write(content)
        

outgroup_gene_names = {'K00370':['dms','tor'],
                       'K00371':['dms','tor'],
                       'K10535':['nrfA','_ONR'],
                       'K10944':['bmo'],
                       'K10945':['bmo'],
                       'K10946':['bmo']}

ref_file = '/home-user/thliao/project/nitrogen_cycle/nitrification/reference_genomes/outgroup and reference.xlsx'
ref_df = pd.read_excel(ref_file,index_col=None)
ref_df = ref_df.loc[ref_df.loc[:,'note']!='removed',:] 
sub_ref_df = ref_df.loc[ref_df.loc[:,'outgroup/ref for which KO']==ko,:]
sub_ref_df = sub_ref_df.loc[sub_ref_df.loc[:,'phylum/class']!='Thaumarchaeota',:]

# step1 write unqiue seqs
fa_file = join(odir,'unique_seqs.faa')
num_unique_seq = len(protein_df.loc[:,'seq'].unique())
sub_protein_df = protein_df.drop_duplicates('seq')
with open(fa_file,'w') as f1:
    for aid,row in sub_protein_df.iterrows():
        seq = row['seq']
        print(f'>{aid}\n{seq}',file=f1)
# filter/ 
# step2 filter false postitive 
ofile = join(odir,f'{ko}.kofam.out')
infa = fa_file[::]
run(f"{kofam_scan} -o {ofile} --cpu 64 -f mapper-one-line --no-report-unannotated {infa} -p /home-user/thliao/data/kofam/profiles/{ko}.hmm")
confirmed_id = [_.strip().split('\t')[0] for _ in open(ofile,'r').read().split('\n') if _]
remained_records = [_ for _ in SeqIO.parse(infa,format='fasta') if _.id in confirmed_id]
new_fa_file = join(odir,'unique_seqs_filtered.faa')
if exists(filter_id_txt):
    ids = [_.strip() for _ in open(filter_id_txt).read().split('\n')]
    remained_records = [_ for _ in remained_records if (_.id not in ids) and (_.id.replace('_',' ') not in ids)]
with open(new_fa_file,'w') as f1:
    SeqIO.write(remained_records,f1,format='fasta-2line')

    
# step3 summarize the distribution of sequences
# lengths distribution
len_seqs = [len(_.seq) for _ in remained_records]
print('25 percentile has %s AA' % pd.np.percentile(len_seqs,25))
print('75 percentile has %s AA' % pd.np.percentile(len_seqs,75))
_s = sorted(remained_records,key=lambda x:len(x.seq))
print('longest seq is %s, has %s AA' % (_s[-1].id,len(_s[-1].seq)))
print('shortest seq is %s, has %s AA' % (_s[0].id,len(_s[0].seq)))
##

# step4 annotate added reference and output and add reference and outgroup into seq
used_records = list(SeqIO.parse(new_fa_file,format='fasta'))
final_records,ref_id2info = add_ref_seq(sub_ref_df,used_records)
ref_id2info,ref_info2style = get_colors_general(ref_id2info)
infa = join(odir, ko+'.fa')
with open(infa,'w') as f1:
    SeqIO.write(final_records,f1,format='fasta-2line')
    
# step5 alignment and build tree
ofile = join(odir, ko+'.aln')
if not exists(ofile):
    print(f'mafft --maxiterate 1000 --genafpair --thread -1 {infa} > {ofile}')#, shell=1)

if not exists( ofile.replace('.aln','.treefile')):
    #pass
    if tree_exe == 'iqtree':
        print(f'iqtree -nt 50 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre {odir}/{ko} -s {ofile}')#,shell=1)
    else:
        n_file = ofile.replace('.aln','.treefile')
        print(f'FastTree {ofile} > {n_file}')#,shell=1)

t = root_tree_with(ofile.replace('.aln','.treefile'),
                    gene_names=outgroup_gene_names.get(ko,[]),
                    format=0)
renamed_tree(t,outfile=ofile.replace('.aln','.sorted.newick'),
                ascending=True)

# generateing annotation files

edl = EntrezDownloader(
    # An email address. You might get blocked by the NCBI without specifying one.
    email='l0404th@gmail.com',
    # An API key. You can obtain one by creating an NCBI account. Speeds things up.
    api_key='ccf9847611deebe1446b9814a356f14cde08',
    num_threads=30,                       # The number of parallel requests to make
    # The number of IDs to fetch per request
    batch_size=500,
    pbar=True                             # Enables a progress bar, requires tqdm package
)
    
remained_records_ids = [_.id for _ in remained_records]
general_df = pd.read_csv(join(odir,'protein2INFO.tab'),sep='\t',index_col=0,low_memory=False)
sub_df = general_df.reindex(remained_records_ids)
# biosample_df = pd.read_excel(join(odir,'biosample2info.xlsx'),index_col=0)
# bioproject_df = pd.read_excel(join(odir,'bioproject2info.xlsx'),index_col=0)
# biosample_df = biosample_df.drop_duplicates().reindex(sub_df.loc[:,'BioSample'])
# bioproject_df = bioproject_df.drop_duplicates().reindex(sub_df.loc[:,'BioProject'])


results, failed = edl.esummary(db='protein',
                                ids=remained_records_ids,
                                result_func=lambda x: Entrez.read(
                                    io.StringIO(x)))
id2tax = {}
id2org = {}
for r in results:
    aid = r['AccessionVersion']
    tid = r['TaxId'].real
    lineage = ncbi.get_lineage(tid)
    rank = ncbi.get_rank(lineage)
    rank = {v: k for k, v in rank.items()}
    names = ncbi.get_taxid_translator(lineage)
    if names.get(rank.get('phylum',''),'ENV') == 'Proteobacteria':
        id2tax[aid] = names.get(rank.get('class',''),'ENV')
    else:
        id2tax[aid] = names.get(rank.get('phylum',''),'ENV')
    id2org[aid] = names[tid]
        
        
        
id2info,info2col = get_colors_general(id2tax,now_info2style= ref_info2style)
id2info.update(ref_id2info)
info2col.update(ref_info2style)
new_text = to_node_symbol(ofile.replace('.aln','.sorted.newick'))
with open(join(odir,f'{ko}_node_bootstrap.txt'),'w') as f1:
    f1.write(new_text)
    

write2colorbranch_clade(id2info,
                        odir,
                        info2col,
                        treefile=ofile.replace('.aln','.sorted.newick'),
                        unique_id=ko,
                        info_name='branch_color',
                        no_legend=False)
        
full_text = to_label(id2org)
with open(join(odir, f'{ko}_label.txt'), 'w') as f1:
    f1.write(full_text)
        