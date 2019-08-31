import pandas as pd
import gffutils
from os.path import join, dirname, basename, exists, abspath
import os
from collections import defaultdict, Counter
from tqdm import tqdm
from scipy.spatial.distance import pdist, squareform


def read_gff(gff_file):
    os.makedirs('./tmp/', exist_ok=True)
    odb_file = join('./tmp/', basename(gff_file) + '.db')
    if not exists(odb_file):
        db = gffutils.create_db(gff_file,
                                odb_file,
                                force=True, )
    else:
        db = gffutils.FeatureDB(odb_file, keep_order=True)
    return db


def get_all_gene_pos(genome_file, CDS_names=None):
    if not exists(genome_file):
        tqdm.write('wrong file path... pass it')
        return None, None
    db = read_gff(genome_file)
    if CDS_names is None:
        all_CDS = [_.id for _ in db.all_features()]
    else:
        all_CDS = CDS_names
    gene2pos = defaultdict(dict)
    order_list = []
    last_end = None
    for fea_id in all_CDS:
        fea = db[fea_id]
        if 'gene' not in fea.id:
            gene2pos[fea_id]['contig_name'] = fea.seqid
            gene2pos[fea_id]['strand'] = fea.strand
            gene2pos[fea_id]['start'] = fea.start
            gene2pos[fea_id]['end'] = fea.end
            gene2pos[fea_id]['previous end'] = last_end
            last_end = fea.end
            order_list.append(fea_id)
    order_list = tuple(order_list)
    return gene2pos, order_list


def get_locus2group(df):
    locus2group = {}
    for group, row in df.iterrows():
        for locus in [locus for _ in row.values
                      for locus in str(_).split(',')]:
            locus2group[locus.split('|')[-1]] = group
    return locus2group


def get_neighbour(target_locus,
                  _order_tuple,
                  locus2group,
                  num_neighbour=20):
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


def determine_locus(row, genome2order_tuple, locus2group):
    locus_all = [(locus.split('|')[-1], genome)
                 for genome, _ in row.items()
                 for locus in str(_).split(',')]
    collect_df = []
    for target_locus, genome in locus_all:
        if target_locus == 'nan':
            continue
        _order_tuple = genome2order_tuple[genome]
        left_n, right_n = get_neighbour(target_locus,
                                        _order_tuple,
                                        locus2group,
                                        num_neighbour=5)
        row_df = pd.DataFrame.from_dict({target_locus: Counter(left_n + right_n)}, orient='index')
        collect_df.append(row_df)
    locus_neighbour_df = pd.concat(collect_df, axis=0)
    return locus_neighbour_df


def group_out():
    pass


def main(infile, prokka_o):
    if not exists(abspath(prokka_o)):
        raise Exception("wrong prokka output directory")
    total_df = pd.read_csv(infile, sep='\t', index_col=0)
    # get all index which contains duplicated genes
    sub_idx = total_df.index[total_df.applymap(lambda x: ',' in str(x)).any(1)]
    tqdm.write('detect %s of duplicated row' % len(sub_idx))
    genomes_files = [join(prokka_o, _, '%s.gff' % _)
                     for _ in total_df.columns]
    genome2gene_info = {}
    genome2order_tuple = {}
    tqdm.write('iterating all gff for collecting positional informations')
    for genome_file in tqdm(genomes_files):
        genome_name = basename(dirname(genome_file))
        _gene_info, _order_tuple = get_all_gene_pos(genome_file)
        if _gene_info is not None:
            genome2gene_info[genome_name] = _gene_info
            genome2order_tuple[genome_name] = _order_tuple
    total_df = total_df.loc[:, list(genome2order_tuple.keys())]
    locus2group = get_locus2group(total_df)
    tqdm.write('start to split duplication...')

    for group_id in sub_idx:
        row = total_df.loc[group_id, :]
        locus_neighbour_df = determine_locus(row, genome2order_tuple, locus2group)



def cli():
    pass
