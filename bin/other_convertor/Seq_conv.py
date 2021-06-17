from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from os.path import *
import os


def gbk2faa(in_gbk,out_faa):
    records = SeqIO.parse(in_gbk,format='genbank')
    
    faa_list = []
    for contig in records:
        for fea in contig.features:
            if fea.type == 'CDS':
                seq = SeqRecord(seq=Seq(fea.qualifiers['translation'][0]),
                            id=fea.qualifiers['locus_tag'][0],
                            name=fea.qualifiers.get('gene',[''])[0],
                            description=fea.qualifiers['product'][0])
                faa_list.append(seq)
    if not exists(dirname(out_faa)):
        os.mkdir(dirname(out_faa))
    with open(out_faa,'w') as f1:
        SeqIO.write(faa_list,f1,'fasta-2line')

# for gbk in glob('./*.gbff'):
#     name = gbk.split('/')[-1].split('_ASM')[0]
#     SeqIO.convert(gbk,in_format='genbank',out_file=f'./{name}.fna'  ,out_format='fasta-2line')


# for _ in tqdm(genomes):
#     in_gbk = f"./gene_arrangments/split_gbk/{_}.gbk"
#     if exists(in_gbk):
#         out_faa = f"./gene_arrangments/split_gbk/faa/{_}.faa"
#         gbk2faa(in_gbk,out_faa)