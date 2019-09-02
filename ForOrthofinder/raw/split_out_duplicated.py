import pandas as pd
import gffutils
from os.path import join, dirname, basename, exists, abspath
import os
from collections import defaultdict, Counter
from tqdm import tqdm
from scipy.spatial.distance import pdist, squareform
from sklearn.neighbors import NearestNeighbors
import click
import pickle
from glob import glob


def preprocess_locus_name(locus):
    locus = str(locus).split('|')[-1]
    return locus


def read_gff(gff_file):
    os.makedirs('./tmp/', exist_ok=True)
    odb_file = join('./tmp/', basename(gff_file) + '.db')
    if not exists(odb_file):
        db = gffutils.create_db(gff_file,
                                odb_file,
                                force=True, sort_attribute_values=True)
    else:
        db = gffutils.FeatureDB(odb_file, keep_order=True, sort_attribute_values=True)
    return db


def get_all_gene_pos(genome_file, CDS_names=None):
    if not exists(genome_file):
        tqdm.write('wrong file path[%s]... pass it' % genome_file)
        return None, None
    db = read_gff(genome_file)
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


def split_out(row, genome2order_tuple, locus2group):
    # core function for splitting the group apart
    # convert neighbours of each target_locus into a counter matrix
    collect_df = []
    for genome, locus in row.items():
        if locus == 'nan' or pd.isna(locus):
            continue
        for locus in locus.split(','):
            target_locus = preprocess_locus_name(locus)
            _order_tuple = genome2order_tuple[genome]
            left_n, right_n = get_neighbour(target_locus,
                                            _order_tuple,
                                            locus2group,
                                            num_neighbour=5)
            if left_n is None or (not left_n and not right_n):
                # None is beacuse it could not find this locus at gff
                # empty is because no neighbour besides it
                continue
            _df = pd.DataFrame.from_dict({target_locus: Counter(left_n + right_n)}, orient='index')
            collect_df.append(_df)
    this_df = pd.concat(collect_df, axis=0, sort=True)
    # calculated the euclidean distance and use nearestNeighbors to get the nearestNeighbors
    # for each locus
    eu_dist = pd.DataFrame(squareform(pdist(this_df.fillna(0))),
                           index=this_df.index,
                           columns=this_df.index)
    model = NearestNeighbors(metric='precomputed', n_neighbors=this_df.shape[0] - 1)
    model.fit(eu_dist.values)
    order_neighbors = model.kneighbors(return_distance=False)
    order_neighbors = pd.np.apply_along_axis(lambda x: this_df.index[x], 0, order_neighbors)
    target_l2neighbours = dict(zip(this_df.index,
                                   order_neighbors))
    # get the name instead of index.
    remained_genomes = set([_.split('_')[0]
                            for _ in target_l2neighbours.keys()])
    group2infos = defaultdict(dict)
    group_num = 1
    while len(remained_genomes) >= 2:
        target_locus = list(target_l2neighbours.keys())[0]
        target_genome = target_locus.split('_')[0]
        target_neighbours = target_l2neighbours.pop(target_locus)
        # use pop, also drop the target_locus
        others_genomes = remained_genomes.difference({target_genome})
        group2infos[group_num][target_genome] = target_locus
        for other_g in others_genomes:
            _cache = [other_l
                      for other_l in target_neighbours
                      if other_l in target_l2neighbours and other_l.startswith(other_g)
                      ]
            # in theory, _cache won't empty?
            group2infos[group_num][other_g] = _cache[0]
            # print(_cache[0]) # debug for
            target_l2neighbours.pop(_cache[0])
            remained_genomes = set([_.split('_')[0]
                                    for _ in target_l2neighbours.keys()])
        group_num += 1

    if len(remained_genomes) == 1:
        genome = list(remained_genomes)[0]
        remained_locus = list(target_l2neighbours)
        for locus in remained_locus:
            group2infos[group_num][genome] = locus
            group_num += 1
    return group2infos


def get_locus2group(df):
    locus2group = {}
    for group, row in df.iterrows():
        for locus in [locus for _ in row.values
                      for locus in str(_).split(',')]:
            locus = preprocess_locus_name(locus)
            locus2group[locus] = group
    return locus2group


def main(infile, prokka_o):
    if not exists(abspath(prokka_o)):
        raise Exception("wrong prokka output directory")
    OG_df = pd.read_csv(infile, sep='\t', index_col=0)
    # get all index which contains duplicated genes
    genomes_files = [glob(join(prokka_o, _ + '*', '%s.gff' % _))[0]
                     for _ in OG_df.columns]
    # use glob and wildcard to capture all real gff files
    if exists('./tmp/genome2gene_info'):
        tqdm.write('detect previous intermediated file, used it to process')
        genome2gene_info = pickle.load(open('./tmp/genome2gene_info', 'rb'))
        genome2order_tuple = pickle.load(open('./tmp/genome2order_tuple', 'rb'))
    else:
        genome2gene_info = {}
        genome2order_tuple = {}
        tqdm.write('iterating all gff for collecting positional information')
        for genome_file in tqdm(genomes_files):
            genome_name = basename(genome_file).replace('.gff', '')
            # use the file name instead of the directory name
            # prokka may remove some string from directory to construct the genome name
            _gene_info, _order_tuple = get_all_gene_pos(genome_file)
            if _gene_info is not None:
                genome2gene_info[genome_name] = _gene_info
                genome2order_tuple[genome_name] = _order_tuple
        # storge temp data
        os.makedirs('./tmp', exist_ok=True)
        pickle.dump(genome2gene_info, open('./tmp/genome2gene_info', 'wb'))
        pickle.dump(genome2order_tuple, open('./tmp/genome2order_tuple', 'wb'))
    OG_df = OG_df.loc[:, list(genome2gene_info.keys())]
    OG_df = OG_df.loc[~OG_df.isna().all(1), :]
    sub_idx = OG_df.index[OG_df.applymap(lambda x: ',' in str(x)).any(1)]
    tqdm.write('detect %s of duplicated row' % len(sub_idx))
    locus2group = get_locus2group(OG_df)
    modify_df = OG_df.copy()
    tqdm.write('collecting all required info, start to split duplicated OG')
    for group_id in tqdm(sub_idx):
        row = OG_df.loc[group_id, :]
        new_group2info = split_out(row, genome2order_tuple, locus2group)
        new_df = pd.DataFrame.from_dict(new_group2info, orient='index')
        new_df.index = [group_id + '_%s' % _
                        for _ in new_df.index]
        new_df = new_df.applymap(lambda x: '%s|%s' % (x.split('_')[0], x) if '_' in str(x) else x)
        modify_df = modify_df.drop(group_id)
        modify_df = modify_df.append(new_df, sort=True)

    # sort the genome with the number of contigs
    order_columns = sorted(modify_df.columns,
                           key=lambda x: len(genome2order_tuple[x]))
    modify_df = modify_df.reindex(columns=order_columns)
    return modify_df


@click.command(help="This script is mainly for splitting duplicated Orthogroup according to the similarity of neighbours.")
@click.option("-i", "infile", help='input file. normally is the concated orthfinder output')
@click.option("-o", "ofile", help='output file')
@click.option("-p", "prokka_dir", help='path of prokka output')
def cli(infile, prokka_dir, ofile):
    modify_df = main(infile, prokka_dir)
    if not dirname(ofile):
        os.makedirs(dirname(ofile))
    modify_df.to_csv(ofile, sep='\t', index=1)


if __name__ == '__main__':
    cli()
