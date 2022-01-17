import io
import multiprocessing as mp
import os
import pickle
from collections import defaultdict, Counter
from glob import glob
from os.path import join, dirname, basename, exists, abspath, expanduser

import click
import gffutils
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.spatial.distance import pdist, squareform
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm

tmp_dir = join(os.environ.get('PWD'), '.tmp')


def process_path(path):
    if not '/' in path:
        path = './' + path
    if path.startswith('~'):
        path = expanduser(path)
    if path.startswith('.'):
        path = abspath(path)
    return path


def preprocess_locus_name(locus):
    locus = locus.strip()
    if '|' in locus:
        locus = str(locus).split('|')[-1].split(' ')[0]
    else:
        locus = str(locus).split(' ')[0]
    locus = locus.strip()
    return locus


def read_gff(gff_file, id_spec):
    os.makedirs(tmp_dir, exist_ok=True)
    odb_file = join(tmp_dir, basename(gff_file) + '.db')
    if not exists(odb_file):
        db = gffutils.create_db(gff_file,
                                odb_file,
                                id_spec=id_spec,
                                force=True,
                                sort_attribute_values=True,
                                merge_strategy='merge')
    else:
        db = gffutils.FeatureDB(odb_file, keep_order=True, sort_attribute_values=True)
    return db


def get_all_gene_pos(genome_file, CDS_names=None, id_spec='ID'):
    """
    major part for reading gff file.
    1. giving a gff file.
    2. iterating all features and retrieve its id. (normally is locus ID)
    3. return gene2pos info like {locus1: {start:1,end:100,strand:'-',previous end:0}}
    3.5. and also record the contig order
    """
    if not exists(genome_file):
        tqdm.write('wrong file path[%s]... pass it' % genome_file)
        return None, None
    if genome_file.endswith('gbk'):
        print("warning, it might use gff to parse gbk file. please pass `-gbk`  ")
    db = read_gff(genome_file, id_spec)
    if CDS_names is None:
        all_CDS = [_.id for _ in db.all_features()]
    else:
        all_CDS = CDS_names
    gene2pos = defaultdict(dict)
    order_contig_list = []
    last_end = None
    contig_list = []
    contig_name = ''
    for fea_id in all_CDS:
        fea = db[fea_id]
        if 'gene' not in fea.id:
            gene2pos[fea_id]['contig_name'] = fea.seqid
            if fea.seqid != contig_name:
                if contig_list:
                    order_contig_list.append(tuple(contig_list))
                contig_list = []
                contig_name = fea.seqid
            gene2pos[fea_id]['strand'] = fea.strand
            gene2pos[fea_id]['start'] = fea.start
            gene2pos[fea_id]['end'] = fea.end
            gene2pos[fea_id]['previous end'] = last_end
            last_end = fea.end
            contig_list.append(fea_id)
    # get final contig
    order_contig_list.append(tuple(contig_list))
    order_contig_list = tuple(order_contig_list)
    return gene2pos, order_contig_list


def process_locus_name_of_gbk(row):
    infos = row.split(' ')
    infos = [_ for _ in infos if _]
    if len(infos) <= 2:
        print(infos)
    try:
        length = infos[1].split('_length_')[1].split('_')[0]
    except:
        import pdb;
        pdb.set_trace()
    name = 'LOCUS' + ' ' * 7 + infos[1].split('_length')[0] + f' {length} bp  ' + '  DNA  linear  20-Jan-2020'
    return name + '\n'


def read_gbk(gbk_file):
    rows = open(gbk_file).readlines()
    rows = [process_locus_name_of_gbk(_)
            if '_length_' in _ and _.startswith('LOCUS') else _
            for _ in rows]
    records = SeqIO.parse(io.StringIO(''.join(rows)), format='genbank')
    records = list(records)
    return records


def get_all_CDS_from_gbk(gbk_file, tag='locus_tag'):
    """
    major part for reading gbk file.
    1. giving a gbk file.
    2. iterating all features and retrieve its id. (normally is locus ID)
    3. return gene2pos info like {locus1: {start:1,end:100,strand:'-',previous end:0}}
    3.5. and also record the contig order
    """
    if not exists(gbk_file):
        tqdm.write('wrong file path[%s]... pass it' % gbk_file)
        return None, None
    contigs = read_gbk(gbk_file)
    order_contig_list = []
    gene2pos = defaultdict(dict)
    last_end = None
    for contig in contigs:
        contig_list = []
        contig_name = contig.id
        all_cds = [fea for fea in contig.features if fea.type == 'CDS']
        for fea in all_cds:
            if tag in fea.qualifiers:
                fea_id = fea.qualifiers[tag][0]
            else:
                fea_id = fea.qualifiers['locus_tag'][0]
            gene2pos[fea_id]['contig_name'] = contig_name
            gene2pos[fea_id]['strand'] = '+' if fea.strand == 1 else '-'
            gene2pos[fea_id]['start'] = int(fea.location.start)
            gene2pos[fea_id]['end'] = int(fea.location.end)
            gene2pos[fea_id]['previous end'] = last_end
            last_end = int(fea.location.end)
            contig_list.append(fea_id)
        order_contig_list.append(tuple(contig_list))
    order_contig_list = tuple(order_contig_list)
    gene2pos = dict(gene2pos)
    return gene2pos, order_contig_list


def get_neighbour(target_locus,
                  _order_tuple,
                  locus2group,
                  num_neighbour=20):
    _order_tuple = [_
                    for _ in _order_tuple
                    if target_locus in _]
    if not _order_tuple:
        return None, None
    _order_tuple = _order_tuple[0]
    target_idx = _order_tuple.index(target_locus)
    l_bound = target_idx - num_neighbour if target_idx >= num_neighbour else 0
    r_bound = target_idx + num_neighbour  # too big would not raise error
    left_n = _order_tuple[l_bound:target_idx]
    right_n = _order_tuple[target_idx + 1:r_bound]
    # get all locus name
    # convert it to group name in order to sort.
    left_n = [locus2group.get(locus, 'not collect')
              for locus in left_n]
    right_n = [locus2group.get(locus, 'not collect')
               for locus in right_n]
    return left_n, right_n


def split_out(row, genome2order_tuple,
              locus2group, remained_bar=True, num_neighbour=5):
    # core function for splitting the group apart
    # convert neighbours of each target_locus into a counter matrix
    # todo: improve the efficiency
    collect_df = []
    locus2genome = {}
    # retrieve neighbour genes and convert in into a OG2number_of_presence table
    row = {genome: locus_raw
           for genome, locus_raw in row.items()
           if not (locus_raw == 'nan' or pd.isna(locus_raw))}
    for genome, locus_raw in row.items():
        for locus in locus_raw.split(','):
            locus = locus.strip()
            target_locus = preprocess_locus_name(locus)
            _order_tuple = genome2order_tuple[genome]
            locus2genome[target_locus] = genome
            left_n, right_n = get_neighbour(target_locus,
                                            _order_tuple,
                                            locus2group,
                                            num_neighbour=num_neighbour)
            if (not left_n and not right_n):
                # None is because it could not find this locus at gff
                # empty is because no neighbour besides it
                _df = pd.DataFrame(index=[target_locus])
            else:
                _df = pd.DataFrame.from_dict({target_locus: Counter(left_n + right_n)}, orient='index')
            collect_df.append(_df)
    this_df = pd.concat(collect_df, axis=0, sort=True)
    # index are the target locus, columns are OG, values represent the number of OG present at the neighbours of target locus
    # calculated the euclidean distance and use nearestNeighbors to get the nearestNeighbors
    # for each locus
    eu_dist = pd.DataFrame(squareform(pdist(this_df.fillna(0))),
                           index=this_df.index,
                           columns=this_df.index)
    model = NearestNeighbors(metric='precomputed',
                             n_neighbors=this_df.shape[0] - 1,
                             n_jobs=-1
                             )
    model.fit(eu_dist.values)
    order_neighbors = model.kneighbors(return_distance=False)
    order_neighbors = np.apply_along_axis(lambda x: this_df.index[x],
                                          0,
                                          order_neighbors)
    target_l2neighbours = dict(zip(this_df.index,
                                   order_neighbors))
    # from target locus to their neighbors ordered with descending distances.
    # get the name instead of index.
    remained_genomes = set([locus2genome[_]
                            for _ in target_l2neighbours.keys()])
    group2infos = defaultdict(dict)
    group_num = 1
    while len(remained_genomes) >= 2:
        target_locus = list(target_l2neighbours.keys())[0]
        target_genome = locus2genome[target_locus]
        target_neighbours = target_l2neighbours.pop(target_locus)
        # use pop, also drop the target_locus
        others_genomes = remained_genomes.difference({target_genome})

        group2infos[group_num][target_genome] = "%s|%s" % (target_genome, target_locus) if remained_bar else target_locus
        for other_g in others_genomes:
            _cache = [other_locus
                      for other_locus in target_neighbours
                      if other_locus in target_l2neighbours and locus2genome[other_locus] == other_g
                      ]
            # in theory, _cache won't empty? so I directly use _cache[0] here...
            group2infos[group_num][other_g] = "%s|%s" % (locus2genome[_cache[0]],
                                                         _cache[0]) if remained_bar else _cache[0]
            target_l2neighbours.pop(_cache[0])
            remained_genomes = set([locus2genome[_]
                                    for _ in target_l2neighbours.keys()])
        group_num += 1

    if len(remained_genomes) == 1:
        genome = list(remained_genomes)[0]
        remained_locus = list(target_l2neighbours)
        for locus in remained_locus:
            group2infos[group_num][genome] = "%s|%s" % (locus2genome[locus], locus) if remained_bar else locus
            group_num += 1
    return group2infos


def get_locus2group(df):
    # specific here
    locus2group = {}
    for group, row in df.iterrows():
        for locus in [locus
                      for _ in row.values
                      for locus in str(_).split(',') if locus]:
            locus = preprocess_locus_name(locus)
            locus2group[locus] = group
    return locus2group


def run(args):
    group_id, genome2order_tuple, row, locus2group, remained_bar, num_neighbour = args
    new_group2info = split_out(row,
                               genome2order_tuple,
                               locus2group,
                               remained_bar=remained_bar,
                               num_neighbour=num_neighbour)
    new_df = pd.DataFrame.from_dict(new_group2info, orient='index')
    new_df.index = [group_id + '_%s' % _
                    for _ in new_df.index]
    return group_id, new_df


def main(infile, prokka_o, use_gbk=False, use_pattern=False, threads=20, num_neighbour=5):
    if not use_pattern:
        if not exists(abspath(prokka_o)):
            raise Exception("wrong prokka output directory")
    OG_df = pd.read_csv(infile, sep='\t', index_col=0, low_memory=False)
    # get all index which contains duplicated genes
    SUFFIX = 'gbk' if use_gbk else 'gff'
    if use_pattern:
        genomes_files = list(glob(prokka_o))
        if len(genomes_files) == 0:
            raise IOError(f"the pattern `{prokka_o}` doesn't exists. Please check it again")
    else:
        genomes_files = [glob(join(prokka_o, _ + '*', f'{_}.{SUFFIX}'))[0]
                         for _ in OG_df.columns
                         if glob(join(prokka_o, _ + '*', f'{_}.{SUFFIX}'))]
    if not genomes_files:
        genomes_files = [join(prokka_o, f'{_}.{SUFFIX}')
                         for _ in OG_df.columns
                         if exists(join(prokka_o, f'{_}.{SUFFIX}'))]
    # use glob and wildcard to capture all real gff files
    if exists(join(tmp_dir, 'genome2gene_info')):
        tqdm.write('detect previous intermediated file, used it to in following analysis')
        genome2gene_info = pickle.load(open(join(tmp_dir, 'genome2gene_info'), 'rb'))
        genome2order_tuple = pickle.load(open(join(tmp_dir, 'genome2order_tuple'), 'rb'))
    else:
        genome2gene_info = {}
        genome2order_tuple = {}
        tqdm.write('iterating all gff for collecting positional information')
        for genome_file in tqdm(genomes_files):
            genome_name = basename(genome_file).rpartition('.')[0]
            if genome_name not in OG_df.columns:
                print(f"{genome_file} not in your table. Just ignore it... Might raise error in the following analysis")
                continue
            # use the file name instead of the directory name
            # prokka may remove some string from directory to construct the genome name
            if use_gbk:
                _gene_info, _order_tuple = get_all_CDS_from_gbk(genome_file)
            else:
                _gene_info, _order_tuple = get_all_gene_pos(genome_file)
            if _gene_info is not None:
                genome2gene_info[genome_name] = _gene_info
                genome2order_tuple[genome_name] = _order_tuple
        # stodge temp data
        os.makedirs(tmp_dir, exist_ok=True)
        pickle.dump(genome2gene_info, open(join(tmp_dir, 'genome2gene_info'), 'wb'))
        pickle.dump(genome2order_tuple, open(join(tmp_dir, 'genome2order_tuple'), 'wb'))

    OG_df = OG_df.loc[:, list(genome2gene_info.keys())]
    OG_df = OG_df.loc[~OG_df.isna().all(1), :]
    # remove the columns which doesn't have provided genomic files.
    if OG_df.applymap(lambda x: '|' in str(x)).any().any():
        # special indicators for personal preference
        remained_bar = True
    else:
        remained_bar = False
    sub_idx = OG_df.index[OG_df.applymap(lambda x: ',' in str(x)).any(1)]
    tqdm.write('detect %s of duplicated row' % len(sub_idx))
    locus2group = get_locus2group(OG_df)
    modify_df = OG_df.copy()
    tqdm.write('collecting all required info, start to split duplicated OG. More duplications would make it slower.')

    # due to the large memory usage of `genoem2order_tuple` and `locus2group`
    # Use manager of multiprocessing to share the large data among all processes
    manager = mp.Manager()
    _d = manager.dict()
    _d2 = manager.dict()
    for k, v in genome2order_tuple.items():
        _d[k] = v
    for k, v in locus2group.items():
        _d2[k] = v
    params = []
    for group_id in sub_idx:
        row = OG_df.loc[group_id, :]
        params.append((group_id,
                       _d,
                       row,
                       _d2,
                       remained_bar,
                       num_neighbour,
                       ))

    tqdm.write('running with multiprocessing ')
    with mp.Pool(processes=threads) as tp:
        r = list(tqdm(tp.imap_unordered(run, params),
                      total=len(params)))
    manager.shutdown()

    modify_df = modify_df.drop([_[0] for _ in r])
    tqdm.write('merge all return results, it may take a while')
    added_df = pd.concat([_[1] for _ in r], axis=0, sort=True)
    added_df = added_df.reindex(columns=modify_df.columns)
    _df = pd.concat([modify_df, added_df], axis=1, sort=True)
    modify_df = modify_df.append(added_df, sort=True)

    # sort the genome with the number of contigs
    order_columns = sorted(modify_df.columns,
                           key=lambda x: len(genome2order_tuple[x]))
    modify_df = modify_df.reindex(columns=order_columns)
    tqdm.write('Done!!! return...')
    return modify_df


@click.command(help="This script is mainly for splitting duplicated Orthogroups according to the similarity of their neighbours.")
@click.option("-i", "infile", help='input file. normally is the orthfinder ouput table')
@click.option("-o", "ofile", help='output file')
@click.option("-p", "prokka_dir", help='path of prokka output')
@click.option("-gbk", "use_gbk", is_flag=True, default=False, help='normally it use gff')
@click.option('-t', "threads", help="number of threads used [20]", default=20)
@click.option("-use_pattern", "use_pattern", is_flag=True, default=False,
              help="If you don't have a direct directory contained prokka output. You could pass a pattern e.g './*/*.gbk' to `-p` to retrieve genomic information. ")
@click.option('-nn', "num_neighbour", default=5, help="number of neighbour genes to consider during splitting duplications. Larger values would slower the program [5]")
def cli(infile, prokka_dir, ofile, use_gbk, use_pattern, threads, num_neighbour):
    modify_df = main(infile, prokka_dir, use_gbk, use_pattern, threads, num_neighbour)
    ofile = process_path(ofile)
    if not dirname(ofile):
        os.makedirs(dirname(ofile))
    modify_df.to_csv(ofile, sep='\t', index=1)
    
# python ~/script/evol_tk/ForOrthofinder/bin/split_out_duplicated.py -i ./data_processing/of_out/Results_Dec25/Orthogroups/Orthogroups.tsv -o ./test.tsv -p ./data_processing/20211224/prokka_o -nn 10

if __name__ == '__main__':
    cli()
