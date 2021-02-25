from dna_features_viewer import GraphicFeature, GraphicRecord
from os.path import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from glob import glob

gbk_indir = "/mnt/maple/thliao/data/NCBI/modified_data/prokka_o"
locus2gene = locus2gene_focal_phyla_finer.copy()
gene2color = {}

## todo: format it
def draw_genome(locus, get_features=False, near_num=5, ignore_locus=False, **kwargs):
    genome = f"GCA_{locus.split('_')[0].replace('v', '.')}"
    gbk_file = join(gbk_indir, genome, genome + '.gbk')
    gene = locus2gene[locus]
    records = list(SeqIO.parse(gbk_file, 'genbank'))
    for contig in records:
        all_cds = [_ for _ in contig.features if _.type in 'CDS']
        num_cds = len(all_cds)
        target_fea_idx = [idx for idx, _ in enumerate(all_cds) if _.qualifiers['locus_tag'][0] == locus]
        if not target_fea_idx:
            continue
        print(num_cds - 1, target_fea_idx)
        l_idx = 0 if target_fea_idx[0] - near_num <= 0 else target_fea_idx[0] - near_num
        r_idx = num_cds if target_fea_idx[0] + near_num > num_cds else target_fea_idx[0] + near_num

        faa_list = []
        for fea in all_cds[l_idx:r_idx]:
            seq = SeqRecord(seq=Seq(fea.qualifiers['translation'][0]),
                            id=fea.qualifiers['locus_tag'][0],
                            name=fea.qualifiers.get('gene', [''])[0],
                            description=fea.qualifiers['product'][0])
            faa_list.append(seq)
        with open(f'/home-user/thliao/project/NOB/gene_arrangments/near_faa/{genome}_{locus}.faa', 'w') as f1:
            SeqIO.write(faa_list, f1, 'fasta-2line')

        return_feas = []
        features = []
        target_strand = 0
        for fea in all_cds[l_idx:r_idx]:
            locus_id = fea.qualifiers['locus_tag'][0]
            if locus_id == locus:
                c = '#9a4848'
                gene = locus2gene[locus_id]
                target_strand = fea.strand

            else:
                c = '#ffd700'
                gene = ''  # fea.qualifiers.get('gene',[''])[0]
                if not gene and locus_id in locusnear2ko:
                    gene = locusnear2ko.get(locus_id, {})
                #                 if not gene:
                #                     og = [og for og,_l in OG2l.items() if locus_id in _l]
                #                     if og:
                #                         gene = og[0]
                if gene in gene2color:
                    c = gene2color[gene]
            return_feas.append((locus_id, gene, fea.strand))

            label = ''
            if gene:
                label = gene
            if not ignore_locus:
                label += '_' + locus_id
            features.append(GraphicFeature(start=fea.location.start.real,
                                           end=fea.location.end.real,
                                           strand=fea.strand,
                                           color=c,
                                           label=label))
        if get_features:
            return return_feas
        record = GraphicRecord(features=features, sequence_length=len(contig.seq)
                               )
        cropped_record = record.crop((record.features[0].start, record.features[-1].end - 1))

        ax, _ = cropped_record.plot(**kwargs)
        if target_strand == -1:
            ax.invert_xaxis()