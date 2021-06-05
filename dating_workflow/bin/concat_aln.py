"""
For concatenating multiple alignment files according to some criteria.

1. name with converted/formatted ID (for those genes annotated from Genbank/refseq genome)
2. input a file designating the grouping information.


"""

import os
import random
from collections import defaultdict
from glob import glob
from os.path import join, exists, dirname, basename

import click
import plotly.graph_objs as go
from Bio import AlignIO, SeqIO
from tqdm import tqdm

from dating_workflow.step_script import process_path, convert_genome_ID, convert_genome_ID_rev


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


def set_partition(f1, name, seq, partition_method):
    if partition_method == 'genes':
        f1.write(f"{name}        {str(seq)}\n")
    elif partition_method == '1,2':
        seq = seq[::3] + seq[1::3]
        f1.write(f"{name}        {str(seq)}\n")
    return str(seq), name


def generate_phy_file(outfile, 
                      record_pos_info, 
                      genome_ids, 
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
        for name, start, end, aln_record in record_pos_info:
            if fill_gaps:
                total_num = len(genome_ids)
            else:
                total_num = len(aln_record)
            if remove_identical:
                _total_num = len(set([str(_.seq) for _ in aln_record]))
                num_identical = total_num - _total_num
                print(f"found {num_identical} identical seq")
                total_num = _total_num
            # total_num = len(aln_record)
            num_seq = len(aln_record)
            length_this_aln = aln_record.get_alignment_length()
            if partition_method == 'genes':
                f1.write(f'{total_num}        {length_this_aln}\n')
            elif partition_method == '1,2':
                length_this_aln -= length_this_aln // 3
                f1.write(f'{total_num}        {length_this_aln}\n')
            used_ids = []
            added_seq = []
            for _ in range(num_seq):
                fid = aln_record[_, :].id
                if name_convertor is not None:
                    fid = name_convertor(fid)
                if str(aln_record[_, :].seq) in set(added_seq) and remove_identical:
                    continue
                if fid in genome_ids:
                    # before _ , should be the converted genome id
                    _seq, _id = set_partition(f1,
                                              name=fid,
                                              seq=aln_record[_, :].seq,
                                              partition_method=partition_method)
                    added_seq.append(str(_seq))
                    used_ids.append(fid)
            if fill_gaps:
                for remained_id in set(genome_ids).difference(set(used_ids)):
                    f1.write(f"{remained_id}        {'-' * length_this_aln}\n")

def get_genomes(genome_list,
                simple_concat=True):

    """
    Accepting a file. 
    It could only contain a column of genome names for simple jobs.
    
    Or
    It could contains multiple lines separated with TAB.
    Besides the first column, the following columns should be the gene names in the alignment.

    Returns:
        dict: name2grouping, maybe empty grouping
    """
    # if genome_list is None:
    #     genome_list = join(indir, 'selected_genomes.txt')
        
    rows = open(genome_list, 'r').read().split('\n')
    
    final_name2grouping = defaultdict(set)
    for row in rows:
        if '\t' not in row:
            name = row if simple_concat else convert_genome_ID(row)
            final_name2grouping[name].add(name)
        else:
            name = row.split('\t')[0]
            final_name2grouping[name].add(row.split('\t')[1])
    return final_name2grouping

def get_genes(indir,suffix,gene_list):
    
    order_seqs = sorted(glob(join(indir, f'*.{suffix}')))
    if gene_list is not None:
        _genes = open(gene_list).read().split('\n')
        if exists(str(gene_list)):
            _genes = open(gene_list).read().split('\n')
            
        elif isinstance(gene_list, str):
            _genes = gene_list.split(',')
        gene_list = [_.strip() for _ in _genes if _]
        order_seqs = [_
                          for _ in order_seqs
                          if basename(_).replace(f'.{suffix}', '') in gene_list]
    return order_seqs

def concat_records(order_seqs,
                   final_name2grouping,
                   g2num_miss=defaultdict(int),
                   suffix='aln',
                   simple_concat=True,
                   ):
    """
    Core function for concatenating.

    Args:
        order_seqs ([type]): ordered alignment files sorted by its name.
        final_name2grouping ([type]): contain final_name to its genes grouping. The grouping maybe empty, and it will auto map the gene name with the final name. 
        suffix ([type]): [description]
        simple_concat ([type]): [description]
        g2num_miss ([type]): only for recording

    Returns:
        [type]: [description]
    """
    record_pos_info = []
    las_pos = 0
    name2record = defaultdict(str)
    tqdm.write('itering all requested files ')
    for idx, aln_file in tqdm(enumerate(order_seqs), 
                              total=len(order_seqs)):
        aln_file_name = basename(aln_file).replace(f'.{suffix}', '')
        aln_record = AlignIO.read(aln_file, format='fasta')
        length_this_aln = aln_record.get_alignment_length()
        # record the partition
        name = "part%s" % int(idx + 1)
        start, end = las_pos + 1, length_this_aln + las_pos
        las_pos = end
        record_pos_info.append((name, start, end, aln_record))
        # done recording
        for final_name,grouping in final_name2grouping.items():
            if not grouping:
                records = [aln_r for aln_r in aln_record
                           if (aln_r.id == final_name and simple_concat) or 
                              (aln_r.id.split('_')[0] == final_name and not simple_concat) ]
            else:
                # designated prefix
                records = [aln_r for aln_r in aln_record
                           if aln_r.id.split('_')[0] in grouping]
            assert len(records) <=1
            if records:
                name2record[final_name] += str(records[0].seq)
            else:
                name2record[final_name] += '-' * length_this_aln
                g2num_miss[aln_file_name] += 1    
    return record_pos_info,name2record


@click.command(help="For concating each aln, if it has some missing part of specific genome, it will use gap(-) to fill it")
@click.option("-i", "indir", help="The input directory which contains all separate aln files")
@click.option("-o", "outfile", default=None, help="Path of output concatenated aln file.")
@click.option("-s", "suffix", default='aln', help="suffix for input files [aln]")
@click.option("-gl", "genome_list", default=None, help="indicate the alternative name or path.")
@click.option("-genel", "gene_list", default=None,
              help="list of gene need to concatenate")
@click.option("-rm_I", "remove_identical", is_flag=True, default=False, 
              help='remove identical sequence for some software like Fasttree. default is not removed')
@click.option("-no_graph", "graph", is_flag=True, default=True, 
              help='generating a graph introducing the number of genes among all genomes. default is generating graph')
@click.option("-no_fill", "fill_gaps", is_flag=True, default=True, help="fill with gaps for genomes doesn't contains this gene. default is filling ")
@click.option("-ct", "concat_type", help='partition or phy or both', default='partition')
@click.option("-p", "partition_method", help='partition with genes or 1st,2nd of codons... (genes|1,2) [genes]', default='genes')
@click.option('-fix_ref', 'fix_refseq', help='fix the name of refseq?', default=False, required=False, is_flag=True)
@click.option('-not_add_prefix', 'not_add_prefix', 
              help='file containing a list of id which do not add prefix as others. ', default=None, required=False)
@click.option('-simple', 'simple_concat', is_flag=True, default=False,
              help='do not perform any name transformation ',  required=False)
def main(indir, 
         outfile, 
         genome_list, 
         gene_list, 
         remove_identical, 
         concat_type, 
         graph, 
         fill_gaps, 
         suffix='aln', 
         fix_refseq=False,
         not_add_prefix=None,
         partition_method='genes',
         simple_concat=False):
    """
    The simple_concat indicate that name in `genome_list` is the genome name.
    If it is False, it indicates that name in `genome_list` is converted/formatted genome name like the prefix of locus.
    """
    if fix_refseq:
        prefix = 'GCF_'
    else:
        prefix = 'GCA_'
    if not_add_prefix is not None:
        not_add_prefix_ids = [_ for _ in open(not_add_prefix).read().split('\n') if _]
    else:
        not_add_prefix_ids = []
    # sampleing the genomes
    final_name2grouping = get_genomes(genome_list,simple_concat)        
    # sampling the gene 
    order_seqs = get_genes(indir,suffix,gene_list)
    
    # init parameters
    g2num_miss = {basename(_).replace(f'.{suffix}', ''): 0 for _ in order_seqs}
    
    
    # concat seqs
    record_pos_info,name2record = concat_records(order_seqs,
                                                 final_name2grouping,
                                                 g2num_miss,
                                                 suffix,
                                                 simple_concat
                                                 )
    print(f"Found {len([k for k,v in g2num_miss.items() if v==0])} backbone genes")
    if outfile is None:
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
                f1.write(f'>{convert_genome_ID_rev(final_name, prefix=prefix,not_add_prefix_ids=not_add_prefix_ids)}\n')
            f1.write(f'{seq}\n')

    if remove_identical:
        remove_identical_seqs(outfile)
    if concat_type.lower() in ['both', 'partition']:
        generate_partition_file(outpartition, record_pos_info)
    if concat_type.lower() in ['both', 'phy']:
        gids = open(genome_list, 'r').read().split('\n')
        if simple_concat:
            name_convertor = lambda x: x
        else:
            name_convertor = lambda x: convert_genome_ID_rev(x,not_add_prefix_ids=not_add_prefix_ids)
        generate_phy_file(outphy, record_pos_info, gids,
                          fill_gaps=fill_gaps,
                          remove_identical=remove_identical,
                          partition_method=partition_method,
                          name_convertor=name_convertor)
    if graph:
        generate_stats_graph(g2num_miss, total=len(gids), ofile=ograph)


if __name__ == '__main__':
    main()
