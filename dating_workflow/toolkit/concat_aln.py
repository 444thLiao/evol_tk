import os

from Bio import AlignIO, SeqIO
from glob import glob
from os.path import join, exists, dirname,basename
import click
import random
from collections import defaultdict
import plotly.graph_objs as go
from tqdm import tqdm

def generate_stats_graph(stats,total,ofile):
    fig = go.Figure()
    ascending_names = sorted(list(stats),
                             key=lambda x:stats[x])
    fig.add_bar(x=list(ascending_names),
                    y=[stats[_] for _ in ascending_names])
                    #text=[f"{stats[_]}/{total}" for _ in ascending_names])
    miss_0_genes = len([k for k,v in stats.items() if v ==0])
    fig.layout.title.text = f'stats graph of all alignment files ({total} genomes, {len(stats)} genes, {miss_0_genes} present at all genomes)'
    fig.write_image(ofile,width=1200,height=900)
    
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

def convert_genome_ID(genome_ID):
    # for GCA_900078535.2
    # it will return
    return genome_ID.split('_')[-1].replace('.','v')

def convert_genome_ID_rev(genome_ID):
    # for 900078535v2
    # it will return
    if not genome_ID.startswith('GCA'):
        if '_' in genome_ID:
            genome_ID = genome_ID.split('_')[0]
        return 'GCA_' + genome_ID.replace('v','.')
    else:
        return genome_ID
    
def generate_partition_file(outfile, record_pos_info):
    with open(outfile, 'w') as f1:
        for name, start, end, _ in record_pos_info:
            f1.write(f"Protein, {name} = {start}-{end} \n")


def generate_phy_file(outfile, record_pos_info, genome_ids,no_fill_gaps=True):
    with open(outfile, 'w') as f1:
        for name, start, end, aln_record in record_pos_info:
            total_num = len(genome_ids)
            # total_num = len(aln_record)
            if no_fill_gaps:
                total_num = len(aln_record)
            num_seq = len(aln_record)
            length_this_aln = aln_record.get_alignment_length()
            f1.write(f'{total_num}        {length_this_aln}\n')
            used_ids = []
            for _ in range(num_seq):
                formatted_gid = aln_record[_, :].id.split('_')[0]
                if formatted_gid in genome_ids:
                    # before _ , should be the converted genome id
                    f1.write(f"{convert_genome_ID_rev(formatted_gid)}        {str(aln_record[_, :].seq)}\n")
                    used_ids.append(formatted_gid)
            if no_fill_gaps:
                continue
            
            for remained_id in set(genome_ids).difference(set(used_ids)):
                f1.write(f"{convert_genome_ID_rev(remained_id)}\n{'-' * length_this_aln}\n")


@click.command(help="For concating each aln, if it has some missing part of specific genome, it will use gap(-) to fill it")
@click.option("-i", "indir", help="The input directory which contains all separate aln files")
@click.option("-o", "outfile", default=None, help="path of outfile. default id in the `-i` directory and named `concat_aln.aln`")
@click.option("-s", "suffix", default='aln')
@click.option("-gl", "genome_list", default=None, help="it will read 'selected_genomes.txt', please prepare the file, or indicate the alternative name or path.")
@click.option("-rm_I", "remove_identical", is_flag=True, default=False)
@click.option("-no_graph", "graph", is_flag=True, default=True)
@click.option("-no_fill", "no_fill_gaps", is_flag=True, default=True)
@click.option("-seed", "seed", help='random seed when removing the identical sequences')
@click.option("-ct", "concat_type", help='partition or phy or both', default='partition')
def main(indir, outfile, genome_list, remove_identical, seed, concat_type, graph,no_fill_gaps,suffix='aln'):
    if genome_list is None:
        genome_list = join(indir, 'selected_genomes.txt')
    with open(genome_list, 'r') as f1:
        gids = f1.read().split('\n')
    gids = [convert_genome_ID(_) for _ in gids]
    # from GCA become locus_tag
    record_pos_info = []
    gid2record = {gid: '' for gid in gids}
    
    las_pos = 0
    order_seqs = sorted(glob(join(indir, f'*.{suffix}')))
    g2num_miss = {basename(_).replace(f'.{suffix}',''):0 for _ in order_seqs}
    tqdm.write('itering all requested files ')
    for idx, aln_file in tqdm(enumerate(order_seqs),total=len(order_seqs)):
        aln_file_name = basename(aln_file).replace(f'.{suffix}','')
        aln_record = AlignIO.read(aln_file, format='fasta')
        length_this_aln = aln_record.get_alignment_length()
        # record the partition
        name = "part%s" % int(idx + 1)
        start, end = las_pos + 1, length_this_aln + las_pos
        las_pos = end
        record_pos_info.append((name, start, end, aln_record))
        # done record
        for gid in gid2record:
            records = [_
                       for _ in aln_record
                       if _.id.split('_')[0] == gid]
            if records:
                gid2record[gid] += str(records[0].seq)
            else:
                gid2record[gid] += '-' * length_this_aln
                
                g2num_miss[aln_file_name] +=1
                
    if outfile is None:
        outfile = join(indir, 'concat_aln.aln')
        outpartition = join(indir, 'concat_aln.partition')
        outphy = join(indir, 'concat_aln.phy')
        ograph = join(indir, 'aln_stats.png')
    else:
        if not '/' in outfile:
            outfile = './' + outfile
        if not exists(dirname(outfile)):
            os.makedirs(dirname(outfile))
        outfile = outfile
        outpartition = join(dirname(outfile), 'concat_aln.partition')
        outphy = join(dirname(outfile), 'concat_aln.phy')
        ograph = join(dirname(outfile), 'aln_stats.png')
        
    with open(outfile, 'w') as f1:
        for gid, seq in gid2record.items():
            f1.write(f'>{convert_genome_ID_rev(gid)}\n')
            f1.write(f'{seq}\n')
            

    if remove_identical:
        remove_identical_seqs(outfile, seed=seed)
    if concat_type.lower() in ['both', 'partition']:
        generate_partition_file(outpartition, record_pos_info)
    if concat_type.lower() in ['both', 'phy']:
        generate_phy_file(outphy, record_pos_info, gids,no_fill_gaps=no_fill_gaps)
    if graph:
        generate_stats_graph(g2num_miss,total=len(gids),ofile=ograph)

if __name__ == '__main__':
    main()
