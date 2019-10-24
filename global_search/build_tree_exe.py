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
from tqdm import tqdm
import sys
ncbi = NCBITaxa()
kofam_scan = '/home-user/thliao/software/kofamscan/exec_annotation'
tree_exe = 'iqtree'

args = sys.argv

infa = args[1]
build_tree_alread = (len(args) >= 3)

# infa = './nr_retrieve_amoC/filtered_by_kegg.faa'
# infa = './nr_retrieve_amoB/filtered_by_kegg.faa'
# infa = './nr_retrieve_amoA/cluster_90'
# infa = './nr_retrieve_nxrB/cluster_95'
# infa = './nr_retrieve_nxrA/cluster_95_filtered_lengths.fa'

gene_info = {'kegg': {'nxrA': 'K00370',
                      'nxrB': 'K00371',
                      'hao': 'K10535',
                      'amoA': 'K10944',
                      'amoB': 'K10945',
                      'amoC': 'K10946'},
             'TIGFAM': {'nxrA': '',
                        'nxrB': '',
                        'hao': 'TIGR01703',
                        'amoA': 'TIGR03080',
                        'amoB': 'TIGR03079',
                        'amoC': 'TIGR03078'}}
# narG: TIGR01580
# narH: TIGR01660
####
odir = dirname(infa)
ko = gene_info['kegg'].get(basename(odir).split('_')[-1], '')
if not ko:
    raise IOError

filter_id_txt = join(odir, 'removed_ids.txt')


def run(cmd):
    check_call(cmd, shell=True)


def add_ref_seq(sub_df, used_records):
    used_ids = [str(_.id) for _ in used_records]
    id2info = {}
    record_need_dropped_ids = []
    for _, row in sub_df.iterrows():
        aa_id = row['AA accession'].strip()
        gene_name = row['gene name']
        seq = row['seq']
        info = row['phylum/class']
        if aa_id in used_ids:
            record_need_dropped_ids.append(aa_id)

        used_records.append(SeqIO.read(io.StringIO(
            f'>{aa_id}_{gene_name}\n{seq}'), format='fasta'))
        id2info[f'{aa_id}_{gene_name}'] = info
    final_records = [_ for _ in used_records if str(
        _.id) not in record_need_dropped_ids]
    return final_records, id2info


color_scheme = {'type': {'NOB': '#e41a1c', 'comammox': '#edc31d',
                         'AOB': '#bad5b9', 'AOA': '#358f0f'},
                'phylum/class': {'Thaumarchaeota': '#358f0f',
                                 'Nitrospirae': '#edc31d',
                                 'Gammaproteobacteria': '#78fce0',
                                 'Chloroflexi': '#e41a1c',
                                 'Betaproteobacteria': '#956cb4',
                                 'Alphaproteobacteria': '#8c613c',
                                 'Actinobacteria': '#11FF11',
                                 'Planctomycetes': '#FF66bb',
                                 }}


def get_colors_general(ID2infos, now_info2style={}):
    _c = color_scheme.copy()
    _c = _c['phylum/class']
    colors = px.colors.qualitative.Dark24 + px.colors.qualitative.Light24
    remained_colors = [c for c in colors if c not in now_info2style.values()]
    info2style = {}
    for v in set(ID2infos.values()):
        if v in _c:
            info2style.update({v: _c[v]})
        else:
            one_color = remained_colors.pop(0)
            info2style.update({v: one_color})
    return ID2infos, info2style


def write2colorbranch_clade(id2info, odir, info2color, treefile, unique_id, info_name='type',
                            **kwargs):
    content = to_color_Clade(
        id2info, info2color, treefile, info_name, **kwargs)
    info_name = info_name.replace('/', '_')
    with open(join(odir, f'{unique_id}_{info_name}_colorbranch_clade.txt'), 'w') as f1:
        f1.write(content)


outgroup_gene_names = {'K00370': ['dms', 'tor'],
                       'K00371': ['dms', ],
                       'K10535': ['nrfA', '_ONR'],
                       'K10944': ['bmo'],
                       'K10945': ['bmo'],
                       'K10946': ['bmo']}

ref_file = '/home-user/thliao/project/nitrogen_cycle/nitrification/reference_genomes/outgroup and reference.xlsx'
ref_df = pd.read_excel(ref_file, index_col=None)
ref_df = ref_df.loc[ref_df.loc[:, 'note'] != 'removed', :]
sub_ref_df = ref_df.loc[ref_df.loc[:, 'outgroup/ref for which KO'] == ko, :]
sub_ref_df = sub_ref_df.loc[sub_ref_df.loc[:,
                                           'phylum/class'] != 'Thaumarchaeota', :]


# filter
remained_records = [_ for _ in tqdm(SeqIO.parse(infa, format='fasta'))]
new_fa_file = join(odir, 'prepared.faa')

with open(new_fa_file, 'w') as f1:
    SeqIO.write(remained_records, f1, format='fasta-2line')


# step3 summarize the distribution of sequences
# lengths distribution
len_seqs = [len(_.seq) for _ in remained_records]
print('25 percentile has %s AA' % pd.np.percentile(len_seqs, 25))
print('75 percentile has %s AA' % pd.np.percentile(len_seqs, 75))
_s = sorted(remained_records, key=lambda x: len(x.seq))
print('longest seq is %s, has %s AA' % (_s[-1].id, len(_s[-1].seq)))
print('shortest seq is %s, has %s AA' % (_s[0].id, len(_s[0].seq)))
##

# step4 annotate added reference and output and add reference and outgroup into seq
used_records = list(SeqIO.parse(new_fa_file, format='fasta'))
final_records, ref_id2info = add_ref_seq(sub_ref_df, used_records)
ref_id2info, ref_info2style = get_colors_general(ref_id2info)
prepared_infa = join(odir, ko+'.fa')
if exists(filter_id_txt):
    ids = [_.strip() for _ in open(filter_id_txt).read().split('\n')]
    final_records = [_
                     for _ in final_records
                     if (_.id not in ids) and (_.id.replace('_', ' ') not in ids)]
with open(prepared_infa, 'w') as f1:
    SeqIO.write(final_records, f1, format='fasta-2line')
print('final prepared fa contains ', len(final_records), ' seqs')

# if not exists(ofile):
used_fa_basename = basename(infa).strip('.') + '_aln.dir'
os.makedirs(join(odir, used_fa_basename), exist_ok=True)
# step5 alignment and build tree
ofile = join(odir, used_fa_basename, ko+'.aln')
# , shell=1)
print(
    f'mafft --maxiterate 1000 --genafpair --thread -1 {prepared_infa} > {ofile}')
# , shell=1)
print(
    f"trimal -in {ofile} -out {ofile.replace('.aln','.trimal')} -automated1 -resoverlap 0.55 -seqoverlap 60")
# if not exists( ofile.replace('.aln','.treefile')):
# pass
if tree_exe == 'iqtree':
    # ,shell=1)
    print(
        f"iqtree -nt 50 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre {ofile.replace('.aln','.iqtree')} -s {ofile.replace('.aln','.trimal')}")
else:
    n_file = ofile.replace('.aln', '.treefile')
    # ,shell=1)
    print(f"FastTree {ofile.replace('.aln','.trimal')} > {n_file}")


if build_tree_alread:
    suffix = args[2]
    os.makedirs(join(odir,
                     used_fa_basename,
                     suffix.strip('.')), exist_ok=True)

    final_suffix = '.sorted.newick'
    t = root_tree_with(ofile.replace('.aln', suffix),
                       gene_names=outgroup_gene_names.get(ko, []),
                       format=0)
    final_ids = list(t.get_leaf_names())
    final_tree = join(odir,
                      used_fa_basename,
                      suffix.strip('.'),
                      basename(ofile).replace('.aln', final_suffix))
    renamed_tree(t,
                 outfile=final_tree,
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
    # general_df = pd.read_csv(join(odir,'protein2INFO.tab'),sep='\t',index_col=0,low_memory=False)
    # sub_df = general_df.reindex(remained_records_ids)
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
        if names.get(rank.get('phylum', ''), 'ENV') == 'Proteobacteria':
            id2tax[aid] = names.get(rank.get('class', ''), 'ENV')
        else:
            id2tax[aid] = names.get(rank.get('phylum', ''), 'ENV')
        id2org[aid] = names[tid]
    id2tax = {k: v for k, v in id2tax.items() if k in final_ids}
    id2org = {k: v for k, v in id2org.items() if k in final_ids}

    id2info, info2col = get_colors_general(
        id2tax, now_info2style=ref_info2style)
    ref_id2info = {k: v for k, v in ref_id2info.items() if k in final_ids}
    id2info.update(ref_id2info)
    info2col.update(ref_info2style)
    new_text = to_node_symbol(final_tree)
    with open(join(dirname(final_tree),
                   f'{ko}_node_bootstrap.txt'), 'w') as f1:
        f1.write(new_text)

    write2colorbranch_clade(id2info,
                            dirname(final_tree),,
                            info2col,
                            treefile=final_tree,
                            unique_id=ko,
                            info_name='branch_color',
                            no_legend=False)

    full_text = to_label(id2org)
    with open(join(dirname(final_tree),
                   f'{ko}_label.txt'), 'w') as f1:
        f1.write(full_text)
