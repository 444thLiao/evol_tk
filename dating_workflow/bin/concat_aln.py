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
            num_seq = len(aln_record)
            length_this_aln = aln_record.get_alignment_length()
            
            # setting the number of sequences
            if fill_gaps:
                # fill unused genomes with gaps, thus the number of sequences should equal to the number of genome_ids
                total_num = len(genome_ids)
            else:
                total_num = len(aln_record)
            # processing the seqs                
            if remove_identical:
                _total_num = len(set([str(_.seq) for _ in aln_record]))
                num_identical = total_num - _total_num
                print(f"found {num_identical} identical seq")
                total_num = _total_num
            
            remaining_seq = []
            for _idx in range(num_seq):
                record_id = aln_record[_idx, :].id
                if name_convertor is not None:
                    record_id = name_convertor(record_id)
                    if not record_id:
                        continue
                    remaining_seq.append(_idx)

            
            texts = []
            added_seq = []
            used_ids = []
            for _idx in range(num_seq):
                record_id = aln_record[_idx, :].id
                sequence = aln_record[_idx, :].seq
                
                if name_convertor is not None:
                    record_id = name_convertor(record_id)
                    if not record_id:
                        # if the record has not been converted into a validate name, remove it
                        continue
                if (str(aln_record[_idx, :].seq) in set(added_seq)) and remove_identical:
                    # if the record has been removed, continue
                    continue
                if record_id in genome_ids:
                    # if the 
                    text = set_partition(name=record_id, # might be the changed one
                                              seq=sequence,
                                              partition_method=partition_method)
                    used_ids.append(text.split(' ')[0])
                    added_seq.append(text.split(' ')[-1].strip('\n'))
                    texts.append(text)
            #   fill unused genomes with gaps                  
            if fill_gaps:
                for remained_id in set(genome_ids).difference(set(used_ids)):
                    # missing genomes
                    texts.append(f"{remained_id}        {'-' * length_this_aln}\n")

            # write out the length
            if partition_method == 'genes':
                f1.write(f'{total_num}        {length_this_aln}\n')
            elif partition_method == '1,2':
                length_this_aln -= length_this_aln // 3
                f1.write(f'{total_num}        {length_this_aln}\n') 
            f1.write(''.join(texts))


def get_genes(indir,suffix,gene_list):
    order_seqs = []
    if ',' in indir:
        for _ in [p.strip() for p in indir.split(',')]:
            seqs = get_genes(_,suffix,gene_list)
            order_seqs.extend(seqs)
        return order_seqs

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
                   name2prefix,
                   g2num_miss=defaultdict(int),
                   suffix='aln',
                   simple_concat=True,
                   ):
    """
    Core function for concatenating.

    Args:
        order_seqs ([type]): ordered alignment files sorted by its name.
        name2prefix ([type]): contain final_name to the prefix of its locus
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
        name = f"{aln_file_name}"
        start, end = las_pos + 1, length_this_aln + las_pos
        las_pos = end
        record_pos_info.append((name, start, end, aln_record))
        # done recording
        for name,prefix in name2prefix.items():
            records = [aln_r 
                        for aln_r in aln_record
                        if (aln_r.id.split('_')[0] in prefix) or (aln_r.id in prefix)]
            if len(records) >1:
                exit("Duplicated IDs are found in the input file")
            if records:
                name2record[name] += str(records[0].seq)
            else:
                name2record[name] += '-' * length_this_aln
                g2num_miss[aln_file_name] += 1    
    return record_pos_info,name2record


@click.command(help="For concating each aln, if it has some missing part of specific genome, it will use gap(-) to fill it")
@click.option("-i", "indir", help="The input directory which contains all separate aln files; Multiple paths could be provieded and separated by comma. ")
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
    
    
    # concat seqs
    record_pos_info,name2record = concat_records(order_seqs,
                                                 name2prefix,
                                                 g2num_miss,
                                                 suffix,
                                                 simple_concat)
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
