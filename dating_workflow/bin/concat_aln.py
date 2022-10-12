"""
For concatenating multiple alignment files according to some criteria.

1. name with converted/formatted ID (for those genes annotated from Genbank/refseq genome)
2. input a file designating the grouping information.


"""

from email.policy import default
import os,io
import random
from collections import defaultdict
from glob import glob
from os.path import join, exists, dirname, basename

import click
import plotly.graph_objs as go
from Bio import AlignIO, SeqIO
from tqdm import tqdm

from dating_workflow.step_script import process_path, convert_genome_ID_rev,get_genomes


def generate_stats_graph(stats, total, ofile):
    fig = go.Figure()
    ascending_names = sorted(list(stats),
                             key=lambda x: stats[x])
    fig.add_bar(x=list(ascending_names),
                y=[stats[_] for _ in ascending_names])
    # text=[f"{stats[_]}/{total}" for _ in ascending_names])
    miss_0_genes = len([k for k, v in stats.items() if v == 0])
    fig.layout.title.text = f'stats graph of all alignment files ({total} genomes, {len(stats)} genes, {miss_0_genes} present at all genomes)'
    fig.write_image(ofile, width=1200, height=900)


def remove_identical_seqs(filename, seed=None):
    if seed is not None:
        random.seed(seed)
    groups_dict = defaultdict(list)
    records = list(SeqIO.parse(filename, format='fasta'))
    ori_length = len(records)
    for record in records:
        groups_dict[str(record.seq)].append(record)
    group_info = open(join(os.path.dirname(filename), "identical_groups.list"), 'w')

    count_ = 0
    new_records = []
    for _, records in groups_dict.items():
        if len(records) == 1:
            new_records.append(records[0])
        else:
            new_records.append(random.choice(records))
            group_name = 'group%s' % count_
            for r in records:
                print("%s\t%s\n" % (group_name, r.id), file=group_info)
            count_ += 1
    if len(new_records) == ori_length:
        print("No identical records found")
    else:
        print("Detect identical records, and now dereplicated into single one. From %s to %s" % (ori_length,
                                                                                                 len(new_records)))
    with open(filename, 'w') as f1:
        SeqIO.write(new_records, f1, format='fasta-2line')


def generate_partition_file(outfile, record_pos_info):
    with open(outfile, 'w') as f1:
        for name, start, end, _ in record_pos_info:
            f1.write(f"Protein, {name} = {start}-{end} \n")

def set_partition(name, seq, partition_method):
    if partition_method == 'genes':
        return f"{name}        {str(seq)}\n"
    elif partition_method == '1,2':
        seq = seq[::3] + seq[1::3]
        return f"{name}        {str(seq)}\n"
    #return str(seq), name


def generate_phy_file(outfile,
                      record_pos_info,
                      genome_ids=[],
                      fill_gaps=True,
                      remove_identical=False,
                      partition_method='genes',
                      name_convertor=None):
    """
    :param outfile:
    :param record_pos_info:
    :param genome_ids: should be the same format id as the record_pos_info, the transforming process should not occur there.
    :param fill_gaps:
    :param remove_identical:
    :param partition_method:
    :return:
    """
    with open(outfile, 'w') as f1:
        if not genome_ids:
            all_genomes = list(set([_.id for ar in record_pos_info for _ in ar[-1] ]))
            genome_ids = all_genomes
            
        for a, b, c, aln_record in record_pos_info:
            num_seq = len(aln_record)
            length_this_aln = aln_record.get_alignment_length()
            total_num = num_seq
            if fill_gaps and num_seq!=len(genome_ids):
                total_num = len(genome_ids)
                
            # remove sequence which are completely same           
            if remove_identical:
                seq2r = {str(s.seq):s for s in aln_record}
                all_seqs = list(seq2r.values())
            else:
                all_seqs = list(aln_record)
            if num_seq - len(all_seqs)!=0:
                print(f"found {num_seq - len(all_seqs)} identical seq")
            final_text = []
            for seq in all_seqs:
                record_id = seq.id
                sequence = str(seq.seq)
                if name_convertor is not None:
                    record_id = name_convertor(record_id)
                    if not record_id: continue
                        # if the record has not been converted into a validate name, remove it
                if record_id in genome_ids:
                    text = set_partition(name=record_id, # might be the changed one
                                         seq=sequence,
                                         partition_method=partition_method)
                    final_text.append(text)
            included_ids = [_.split(' ')[0] for _ in final_text]
            #   fill unused genomes with gaps
            if fill_gaps:
                for remained_id in set(genome_ids).difference(set(included_ids)):
                    # missing genomes
                    final_text.append(f"{remained_id}        {'-' * length_this_aln}\n")
            # write out the length
            if partition_method == 'genes':
                f1.write(f'{total_num}        {length_this_aln}\n')
            elif partition_method == '1,2':
                length_this_aln -= length_this_aln // 3
                f1.write(f'{total_num}        {length_this_aln}\n')
            f1.write(''.join(final_text))


def get_genes(indir,suffix,gene_list=None):
    """
    return a dict with Gene to partition. If no paritions are assigned, each gene to each gene.
    """
    if gene_list is None:
        order_seqs = {}
        if ',' in indir:
            # if multiple directories passed.
            for _ in [p.strip() for p in indir.split(',')]:
                seqs = get_genes(_,suffix,gene_list)
                order_seqs.update(seqs)
            return order_seqs
        ## if not gene_list, return and end
    order_seqs = {fname:basename(fname).replace(f'.{suffix}', '') for fname in sorted(glob(join(indir, f'*.{suffix}')))}
    if gene_list is not None:
        g2p = {}
        if exists(str(gene_list)):
            for row in open(gene_list).read().strip().split('\n'):
                if '\t' in row:
                    g2p[row.split('\t')[0]] = row.split('\t')[1]
                else:
                    g2p[row.split('\t')[0]] = row.split('\t')[1]

        elif isinstance(gene_list, str) and not exists(str(gene_list)):
            g2p.update({k:k for k in gene_list.split(',')})

        order_seqs = {fname: g2p.get(name,
                                     name)
                      for fname,name in order_seqs.items()
                      if name in list(g2p)}
        # if there are second colum in the gene_list file, it will assign a partition name.
    return order_seqs

def concat_records(order_seqs,
                   name2prefix,
                   g2num_miss=defaultdict(int)
                   ):
    """
    Core function for concatenating.

    Args:
        order_seqs ([type]): ordered alignment files sorted by its name.
        name2prefix ([type]): contain final_name to the prefix of its locus
        suffix ([type]): [description]
        g2num_miss ([type]): only for recording

    Returns:
        1. partition information
        2. single concatenated alignment
    """
    partition2aln_files = defaultdict(list)
    for aln_file,partition in order_seqs.items():
        partition2aln_files[partition].append(aln_file)
    # maybe each parition has a single aln_file
    # if no partition is given, it should be name of the alignment file (without suffix).
    fullname2record = defaultdict(str)
    tqdm.write('itering all requested files ')
    record_pos_info = []
    for partition,aln_files in tqdm(partition2aln_files.items(),total=len(partition2aln_files)):
        name2record = defaultdict(str)
        for aln_file in aln_files:
            aln_record = AlignIO.read(aln_file, format='fasta')
            length_this_aln = aln_record.get_alignment_length()
            for final_name,prefix in name2prefix.items():
                records = [aln_r
                            for aln_r in aln_record
                            if (aln_r.id.split('_')[0] in prefix) or (aln_r.id in prefix)]
                if len(records) >1:
                    exit("Duplicated IDs are found in the input file")
                if records:
                    name2record[final_name] += str(records[0].seq)
                else:
                    name2record[final_name] += '-' * length_this_aln
                    g2num_miss[final_name] += 1
        # init the positional information
        if not record_pos_info:
            start = 1
            end = len(name2record[final_name])
        else:
            start = record_pos_info[-1][2] +1
            end = len(name2record[final_name]) + start
        # init a alignment records
        c = ''
        for k,v in name2record.items():
            c+=f">{k}\n{v}\n"
            fullname2record[k]+=v
        aln_record = AlignIO.read(io.StringIO(c), format='fasta')

        record_pos_info.append((partition, start, end, aln_record))
    return record_pos_info,fullname2record


@click.command(help="For concating each aln, if it has some missing part of specific genome, it will use gap(-) to fill it")
@click.option("-i", "indir", help="The input directory which contains all separate aln files; Multiple paths could be provieded and separated by comma. ")
@click.option("-o", "outfile", default=None, help="Path of output concatenated aln file.")
@click.option("-s", "suffix", default='aln', help="suffix for input files [aln]")
@click.option("-gl", "genome_list", default=None, help="file containing a list of species name. If there are two columns (tab), the second should be the name of gene in the alignment file.")
@click.option("-genel", "gene_list", default=None,
              help="file containing a list of gene name (mostly the name of alignment file). If there are two columns (tab), the second should be the name of partition used for partition. If only a single column, all will be treated as a single partition. ")
@click.option("-rm_I", "remove_identical", is_flag=True, default=False,
              help='remove identical sequence for some software like Fasttree. default is not removed')
@click.option("-no_graph", "graph", is_flag=True, default=False,
              help='generating a graph introducing the number of genes among all genomes. default is generating graph')
@click.option("-no_fill", "fill_gaps", is_flag=True, default=True, help="fill with gaps for genomes doesn't contains this gene. default is filling ")
@click.option("-ct", "concat_type", help='partition or phy or both', default='partition')
@click.option("-p", "partition_method", help='partition with genes or 1st,2nd of codons... (genes|1,2) [genes]. If partition information are specificed in the gene_list, this parameter will be useless. ', default='genes')
@click.option('-fix_ref', 'fix_refseq', help='fix the name of refseq?', default=False, required=False, is_flag=True)
@click.option('-simple', 'simple_concat', is_flag=True, default=False,
              help='do not perform any name transformation ',  required=False)
def main(indir,
         outfile,
         genome_list,
         gene_list,
         concat_type,
         graph,
         fill_gaps,
         suffix='aln',
         fix_refseq=False,
         remove_identical=False,
         partition_method='genes',
         simple_concat=True):
    """
    The simple_concat indicate that name in `genome_list` is the genome name.
    If it is False, it indicates that name in `genome_list` is converted/formatted genome name like the prefix of locus.
    """
    if fix_refseq:
        prefix = 'GCF_'
    else:
        prefix = 'GCA_'
    # sampleing the genomes
    name2prefix = get_genomes(genome_list,simple_concat)
    # sampling the gene
    order_seqs = get_genes(indir,suffix,gene_list)

    # init parameters
    g2num_miss = {basename(_).replace(f'.{suffix}', ''): 0 for _ in order_seqs}

    # concat seqs (according to given partitions or a single partition )
    record_pos_info,name2record = concat_records(order_seqs,
                                                 name2prefix,g2num_miss)
    print(f"Found {len([k for k,v in g2num_miss.items() if v==0])} backbone genes")
    if outfile is None and ',' not in indir:
        outfile = join(indir, 'concat_aln.aln')
        outpartition = join(indir, 'concat_aln.partition')
        outphy = join(indir, 'concat_aln.phy')
        ograph = join(indir, 'aln_stats.png')
    else:
        outfile = process_path(outfile)
        if not exists(dirname(outfile)):
            os.makedirs(dirname(outfile))
        outpartition = outfile.rpartition('.')[0] + '.partition'
        outphy = outfile.rpartition('.')[0] + '.phy'
        ograph = join(dirname(outfile), 'aln_stats.png')

    with open(outfile, 'w') as f1:
        for final_name, seq in name2record.items():
            if set(str(seq)) == {'-'}:
                print(f"{final_name} contains only gaps or missing data ")
                continue
            if simple_concat:
                f1.write(f">{final_name}\n")
            else:
                f1.write(f'>{convert_genome_ID_rev(final_name, prefix=prefix)}\n')
            f1.write(f'{seq}\n')

    if remove_identical:
        remove_identical_seqs(outfile)
    if concat_type.lower() in ['both', 'partition']:
        generate_partition_file(outpartition, record_pos_info)
    if concat_type.lower() in ['both', 'phy']:
        gids = list(name2prefix)
        def name_convertor(x):
            tmp = [k for k,v in name2prefix.items() if x.split('_')[0] in v or x in v]
            if not tmp:
                return
            else:
                return tmp[0]
        generate_phy_file(outphy,
                          record_pos_info,
                          gids,
                          fill_gaps=fill_gaps,
                          remove_identical=remove_identical,
                          partition_method=partition_method,
                          name_convertor=name_convertor)
    if graph:
        generate_stats_graph(g2num_miss, total=len(gids), ofile=ograph)


if __name__ == '__main__':
    main()
