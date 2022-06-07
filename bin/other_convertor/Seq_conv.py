from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from os.path import *
import os


def gbk2faa(in_gbk, out_faa, prefix="", include_tRNA=False,return_DNA=False):
    records = list(SeqIO.parse(in_gbk, format="genbank"))
    faa_list = []
    for contig in records:
        idx = 1
        for fea in contig.features:
            if prefix:
                seq_name = f"{prefix}_" + "{:05d}".format(idx)
            else:
                seq_name = fea.qualifiers.get("locus_tag", [""])[0]
            if include_tRNA:
                type_list = ["CDS", "tRNA"]
            else:
                type_list = ["CDS"]

            if fea.type in type_list:
                if return_DNA:
                    _seq = Seq(str(fea.extract(contig).seq))
                else:
                    _seq = Seq(fea.qualifiers.get("translation",[''])[0])
                new_record = SeqRecord(
                    seq=_seq,
                    id=seq_name,
                    name=fea.qualifiers.get("gene", [""])[0],
                    description=fea.qualifiers["product"][0],
                )
                faa_list.append(new_record)
                idx += 1
    if not exists(dirname(out_faa)):
        os.mkdir(dirname(out_faa))
    with open(out_faa, "w") as f1:
        SeqIO.write(faa_list, f1, "fasta-2line")


# for gbk in glob('./*.gbff'):
#     name = gbk.split('/')[-1].split('_ASM')[0]
#     SeqIO.convert(gbk,in_format='genbank',out_file=f'./{name}.fna'  ,out_format='fasta-2line')


# for _ in tqdm(genomes):
#     in_gbk = f"./gene_arrangments/split_gbk/{_}.gbk"
#     if exists(in_gbk):
#         out_faa = f"./gene_arrangments/split_gbk/faa/{_}.faa"
#         gbk2faa(in_gbk,out_faa)
