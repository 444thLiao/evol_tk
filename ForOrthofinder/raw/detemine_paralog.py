import pandas as pd
from os.path import join, dirname, basename, exists
from Bio import SeqIO
from glob import glob

def get_file_name(name, prokka_dir, suffix='faa'):
    path_file = glob(join(prokka_dir, name+'*', name + '.' + suffix))[0]
    # .1 is not robust
    if exists(path_file):
        return path_file
    else:
        raise Exception(path_file, " doesn't exist. Please check it.")


def getOrder(infile, format='fasta'):
    # For now ,don't consider the differ contig problems
    # because the phage doesn't have different contig
    if format == 'fasta':
        # fixme: doesn't consider different contig problems
        # todo: 
        records = SeqIO.parse(infile, format='fasta')
        return tuple([_.id for _ in records])
    elif format == 'gff':
        pass


def get_neighbors(order_locus, sc, num):
    sc_pos = order_locus.index(sc)
    lefts = order_locus[sc_pos - num:sc_pos]
    rights = order_locus[sc_pos:sc_pos + num]
    return lefts, rights


def calculate_p(lefts1, rights1,
                lefts2, rights2):
    intersec_num_l = len(set(lefts1).intersection(set(lefts2)))
    union_num_l = len(set(lefts1).union(set(lefts2)))
    intersec_num_r = len(set(rights1).intersection(set(rights2)))
    union_num_r = len(set(rights1).union(set(rights2)))

    p_l = intersec_num_l/union_num_l
    p_r = intersec_num_r/union_num_r
    return (p_l+p_r)/2

def choose_paralogs(single_copy, paralogs, g2order_locus,
                    locus2g,
                    support_num=2, support_p=0.5):
    if len(single_copy) == 0:
        # time consumed operation
        pass
    else:
        # confirm this way is ok for single copy
        last_lefts = last_rights = []
        for g, sc in single_copy.items():
            order_locus = g2order_locus.get(g, None)
            if order_locus is not None:
                lefts, rights = get_neighbors(order_locus, sc, support_num)
                lefts, rights = tuple([locus2g.get(_) for _ in lefts]), \
                                tuple([locus2g.get(_) for _ in rights])
                if last_lefts and last_rights:
                    cal_p = calculate_p(lefts, rights,
                                        last_lefts, last_rights)
                    print(cal_p)
                    if cal_p <= support_p:
                        raise Exception("weird.....")
                last_lefts = lefts
                last_rights = rights


def get_locus2g(data):
    locus2g = {}
    for og, row in data.iterrows():
        for each in row:
            if not pd.isna(each):
                locus2g.update({l: og
                                for l in str(each).split(',')})
    return locus2g


def main(infile, prokka_dir):
    data = pd.read_csv(infile, sep='\t', index_col=0, low_memory=False)
    data = data.iloc[:,:-10]
    contain_paralog = data.applymap(lambda x: ',' in str(x))
    OG_with_paralog = data.index[contain_paralog.any(1)]
    # if have COMMA, it is a og with paralog
    g2name = {g: get_file_name(g, prokka_dir, 'faa')
              for g in data.columns}
    g2order_locus = {g: getOrder(g2name[g], format='fasta')
                     for g in data.columns}
    locus2g = get_locus2g(data)
    for OG in OG_with_paralog:
        row = data.loc[OG, :]
        notnan_row = row[~row.isna()]
        single_copy = notnan_row[~notnan_row.str.contains(',')]
        paralogs = notnan_row[notnan_row.str.contains(',')]
        assert len(paralogs) > 0
        result = choose_paralogs(single_copy, paralogs, g2order_locus,
                                 locus2g)

if __name__ == "__main__":
    main("/home-user/thliao/project/cyanophage/within_cyanophage/ortho_test/data/Results_Aug09/WorkingDirectory/Orthogroups_3.csv",
         "/home-user/thliao/project/cyanophage/within_cyanophage/prokka_o")