from tqdm import tqdm
from Bio import SeqIO
from api_tools.IO_for.read import _parse_hmmscan

ko2genes_file = "/mnt/home-db/pub/protein_db/kegg/v2019/ko/ko_genes.list"
ko2genes = {}
nifH_genes = []
vnfH_genes = []
for row in tqdm(open(ko2genes_file)):
    knumber, gene = row.split('\t')
    if "K02588" in knumber:
        nifH_genes.append(gene.strip('\n'))
    if "K22899" in knumber:
        vnfH_genes.append(gene.strip('\n'))

target_faa = "/mnt/home-db/pub/protein_db/kegg/v2019/species_prokaryotes.pep"
records = SeqIO.parse(target_faa,format='fasta')
nifH_genes_faa = []
vnfH_genes_faa = []
for r in tqdm(records,total=13405245):
    if r.id in nifH_genes:
        nifH_genes_faa.append(r)
    if r.id in vnfH_genes:
        vnfH_genes_faa.append(r)

# In [25]: print(len(nifH_genes),len(vnfH_genes))
# 875 27

# In [50]: print(len(nifH_genes_faa),len(vnfH_genes_faa))
# 718 15
# some of the genes actually don't stodge inside the pep file.


nif_tab = "/home-user/jjtao/Rhizobiales/kegg_hmmsearch/annotated_tab/K02588.tab"
vnf_tab = "/home-user/jjtao/Rhizobiales/kegg_hmmsearch/moved_ko/K22899.tab"


nif_ = _parse_hmmscan(nif_tab)
vnf_ = _parse_hmmscan(vnf_tab)

nif_ = {k: sorted(v,key=lambda x: x[-1])[0] if len(v)!=0 else []
        for k,v in nif_.items()}
vnf_ = {k: sorted(v,key=lambda x: x[-1])[0] if len(v)!=0 else []
        for k,v in vnf_.items()}
genomes_list = set(list(nif_) + list(vnf_))

vnf_st_nif_genomes = []
nif_st_vnf_locus = []
vnf_st_nif_locus = []
for genome in genomes_list:
    vnf_g = vnf_.get(genome, [])
    nif_g = nif_.get(genome, [])
    if len(vnf_g) != 0 and len(nif_g) != 0:
        if vnf_g[-1] <= nif_g[-1] and vnf_g[0] == nif_g[0]:
            vnf_st_nif_genomes.append(genome)
            vnf_st_nif_locus.append((vnf_g[-1],vnf_g[0]))
        if vnf_g[-1] >= nif_g[-1] and vnf_g[0] == nif_g[0]:
            nif_st_vnf_locus.append((nif_g[-1],vnf_g[0]))


total_faa = "/home-user/jjtao/Rhizobiales/kegg_hmmsearch/concat346_rmdup.faa"
records = list(SeqIO.parse(total_faa,format='fasta'))
vnf_st_nif_locus_faa = [_ for _ in records if _.id in vnf_st_nif_locus]
nif_st_vnf_locus_faa = [_ for _ in records if _.id in nif_st_vnf_locus]

print(len(vnf_st_nif_locus),len(vnf_st_nif_locus_faa),len(nif_st_vnf_locus_faa))

with open('./vnfH_like_nifH.faa','w') as f1:
    SeqIO.write(vnf_st_nif_locus_faa,f1,format='fasta-2line')

with open('./vnfH_genes.faa','w') as f1:
    SeqIO.write(vnfH_genes_faa,f1,format='fasta-2line')

total_r = nifH_genes_faa+vnfH_genes_faa+vnf_st_nif_locus_faa+nif_st_vnf_locus_faa
total_r_new = []
for _ in total_r:
    _.id = _.id.replace(':','_')
    total_r_new.append(_)
with open("./total.faa","w") as f1:
    SeqIO.write(total_r_new,f1,format="fasta-2line")

print(len(nifH_genes_faa),len(vnfH_genes_faa),len(vnf_st_nif_locus_faa),len(nif_st_vnf_locus_faa))
# mafft --auto --thread -1 ./total.faa > ./total.aln
## trimal -in ./total.aln -out ./total.trimal -automated1 -resoverlap 0.55 -seqoverlap 60 -keepheader
# FastTreeMP ./total.aln > ./total.newick

from api_tools.itol_func import *
id2s = {}
for _ in nifH_genes_faa:
    id2s[_.id.replace(':','_')] = 'nifH'
for _ in vnfH_genes_faa:
    id2s[_.id.replace(':','_')] = 'vnfH'
for _ in vnf_st_nif_locus_faa:
    id2s[_.id.replace(':','_')] = 'vnfH_like_nifH'
for _ in nif_st_vnf_locus_faa:
    id2s[_.id.replace(':', '_')] = 'nifH_like_nifH'

text = to_matrix_shape(id2s,"explore nifH")
with open('./itol.txt','w') as f1:
    f1.write(text)