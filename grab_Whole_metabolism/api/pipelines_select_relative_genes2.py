import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
from Bio import SeqIO
from collections import Counter,defaultdict
import os
from os.path import basename, dirname, abspath, exists, join
from subprocess import check_call


def output_seq_from_df(df, oseq):
    with open(oseq, 'w') as f1:
        for _, row in df.iterrows():
            aa_seq = row['AA seq']
            name = row["locus_name"]
            f1.write(f'>{name}\n')
            f1.write(f"{aa_seq}\n")


def run_cmd(cmd):
    check_call(cmd, shell=True)


def main(infile, target_fa, oseq, project_name):
    df = pd.read_csv(infile, sep='\t')
    odir = dirname(oseq)
    if not exists(odir):
        os.makedirs(odir)
    output_seq_from_df(df, oseq)
    dbname = abspath(oseq).rpartition('.')[0]

    # name
    o1_tab = join(odir, 'first_diamond.out')
    intermedia_faa = join(odir, 'first_extract_seq.faa')
    o2_tab = join(odir, 'second_diamond.out')
    o2_tab = '/home-user/thliao/data/metagenomes/new_grab/output/first_extract_seq_KOfam.out'
    final_faa = join(odir, 'confirmed.faa')

    run_cmd(f'diamond makedb --in {oseq} --db {dbname}')
    run_cmd(f"diamond blastp -q {target_fa} -o {o1_tab} -d {dbname} -p 0 -b 5 -c 2")

    tmp_df = pd.read_csv(f'{o1_tab}', sep='\t', header=None)
    records = SeqIO.parse(f'{target_fa}', format='fasta')
    used_gids = set(tmp_df.iloc[:, 0])
    collcect_records = []
    for record in tqdm(records):
        if record.id in used_gids:
            collcect_records.append(record)
    with open(f'{intermedia_faa}', 'w') as f1:
        SeqIO.write(collcect_records, f1, format='fasta-2line')

    run_cmd(f"/home-user/thliao/software/kofamscan/exec_annotation {intermedia_faa} -o {o2_tab} --cpu 2 -f mapper-one-line --no-report-unannotated")

    pre_df = pd.read_csv(f'{o1_tab}', sep='\t', header=None)
    query_df = df.copy()
    query_df = query_df.drop_duplicates('locus_name')
    query_df = query_df.set_index('locus_name')
    get_df = query_df.reindex(pre_df.loc[:, 1])
    pre_df.loc[:, 'cover ratio'] = pre_df.loc[:, 3].values / get_df.loc[:, 'AA seq'].str.len().values
    pre_df.loc[:, 'KO'] = get_df.loc[:, 'Orthology(single)'].values
    pre_df.loc[:, 'KO name'] = get_df.loc[:, 'KO name'].values
    aft_str = open(f'{o2_tab}', 'r').read().split('\n')
    aft_dict = dict([(row.split('\t')[0],
                      row.split('\t')[1:]) for row in aft_str])
    aft_list_ID = list(aft_dict.keys())

    not_in_df = pre_df.loc[~pre_df.loc[:, 0].isin(aft_list_ID), :]
    in_df = pre_df.loc[pre_df.loc[:, 0].isin(aft_list_ID), :]
    collect_seq = defaultdict(dict)  # confirmed seq
    # for in_df
    for _, row in tqdm(in_df.iterrows()):
        locus = row[0]
        ko_name = row['KO name']
        ko1 = row["KO"]
        ko2 = aft_dict[locus]
        if ko1 not in ko2 and len(ko2)>=1:
            pass
        elif ko1 in ko2 and len(ko2) == 1:
            collect_seq[locus]['real KO'] = ko1
            collect_seq[locus]['KO name'] = ko_name
        elif ko1 in ko2 and ko2.index(ko1) != 0:
            collect_seq[locus]['real KO'] = ko2[0]
            collect_seq[locus]['is paralog'] = 'TRUE'
            collect_seq[locus]['likely KO'] = ko1


    print("contains %s confirmed sequence" % len(real_N_metabolism_genes))
    real_N_metabolism_genes = set(real_N_metabolism_genes)
    records = SeqIO.parse(f'{intermedia_faa}', format='fasta')
    collect_reads = [_ for _ in records if _.id in set(real_N_metabolism_genes)]

    with open(f'{final_faa}', 'w') as f1:
        SeqIO.write(collect_reads, f1, format='fasta-2line')

    locus_list = [_.id for _ in collect_reads]
    locus2gene_df = pd.read_csv(infile, sep='\t', index_col=0)
    sub_df = pre_df.loc[pre_df.loc[:, 0].isin(locus_list), :]

    locus2ko = dict()
    locus2module = dict()
    locus2completeOrthos = dict()

    choose_highest_one = sub_df.sort_values([0, 10]).drop_duplicates(0)
    for rid, row in tqdm(choose_highest_one.iterrows(),
                         total=choose_highest_one.shape[0]):
        locus_tag = row[1]
        seq_name = row[0]
        name = locus2gene_df.loc[locus_tag, "Name"]
        ortho = locus2gene_df.loc[locus_tag, "Orthology(single)"]
        module = locus2gene_df.loc[locus_tag, "module Name"]
        completeOrthos = locus2gene_df.loc[locus_tag, "Orthology(total)"]

        locus2ko[seq_name] = ortho
        locus2module[seq_name] = module
        locus2completeOrthos[seq_name] = completeOrthos

    locus2name = {}
    # module_counts = Counter([tuple(sorted(_)) for _ in locus2module.values()])
    for locus, ko in tqdm(locus2ko.items()):
        if not isinstance(ko, str):
            ko = ';'.join(set(list(ko.values)))
        _sub_df = locus2gene_df.loc[locus2gene_df.loc[:, 'Orthology(single)'] == ko, 'KO name']
        name = ';'.join(list(set(list(_sub_df))))
        locus2name[locus] = name

    sample2info = pd.read_csv('/home-user/thliao/data/metagenomes/concat_all/sample2infos.tsv', sep='\t', header=0, index_col=1)
    locus2info_df = pd.DataFrame(columns=["locus",
                                          "locus_prefix",
                                          'sample name',
                                          'source project',
                                          'ko(single)',
                                          'Gene name(N metabolism)',
                                          'ko(complete)',
                                          'module',
                                          ])
    from collections import defaultdict

    count_ = 0
    locus2info_dict = defaultdict(dict)
    for locus, g_name in tqdm(locus2name.items(),
                              total=len(locus2name)):
        locus_prefix = locus.split('_')[0]
        sname, source = sample2info.loc[locus_prefix, :].values[:2]
        if not isinstance(locus2ko[locus], str):
            for rid, v in enumerate(locus2ko[locus]):
                locus2info_df.loc[count_, :] = [locus,
                                                locus_prefix,
                                                sname,
                                                source,
                                                locus2ko[locus].values[rid],
                                                g_name,
                                                locus2completeOrthos[locus].values[rid],
                                                locus2module[locus].values[rid]]
                count_ += 1
        else:
            locus2info_df.loc[count_, :] = [locus,
                                            locus_prefix,
                                            sname,
                                            source,
                                            locus2ko[locus],
                                            g_name,
                                            locus2completeOrthos[locus],
                                            locus2module[locus]]
            count_ += 1

    order_sample2info = sample2info.copy()
    order_sample2info = order_sample2info.reindex([_.split('_')[0] for _ in locus2info_df.locus])
    c = ['superkingdom',
         'phylum',
         'class',
         'order',
         'family',
         "genus",
         'species',
         'superkingdom(from metadata)',
         'phylum(from metadata)',
         'class(from metadata)',
         'order(from metadata)',
         'family(from metadata)',
         'genus(from metadata)',
         'species(from metadata)', ]
    locus2info_df = locus2info_df.reindex(columns=list(locus2info_df.columns) + c)
    locus2info_df.loc[:, c] = order_sample2info.loc[:, c].values
    locus2info_df.to_csv(join(odir, 'confirmed_locus2info.tsv'), sep='\t', index=0)
    ############################################################
    from collections import defaultdict
    from bioservices import KEGG
    def count_method(df, level):
        each_v2num = {}
        each_v2snames = defaultdict(list)
        for each_v in df.loc[:, level].unique():
            _df = df.loc[df.loc[:, level] == each_v, :]
            num_count = 0
            counted_samples = []
            for ko_complete in _df.loc[:, 'ko(complete)'].unique():
                if ',' not in ko_complete:
                    _cache = _df.loc[_df.loc[:, 'ko(complete)'] == ko_complete, 'sample name'].unique()
                    num_count += len([_ for _ in _cache if _ not in counted_samples])
                    counted_samples += list(_cache)
                    continue

                for sname in _df.loc[:, 'sample name'].unique():
                    if sname in counted_samples:
                        continue
                    sub_df = _df.loc[_df.loc[:, 'sample name'] == sname, :]
                    if not set(ko_complete.split('+')).difference(set(sub_df.loc[:, 'ko(single)'].values)):
                        counted_samples.append(sname)
                        num_count += 1
            each_v2snames[each_v] = counted_samples
            each_v2num[each_v] = num_count
        each_v2num = pd.DataFrame.from_dict({0: each_v2num}).loc[:, 0].sort_values(ascending=False)
        return each_v2num, each_v2snames

    kegg = KEGG()
    ko2info = pd.read_csv('./ko_info.csv', sep='\t', index_col=0)
    ko2gname = dict(zip(ko2info.index,
                        ko2info.loc[:, 'gene name']))
    used_ko = locus2info_df.loc[:, 'ko(single)']
    koSingle2name = {}
    rename_dict = {'K20934': 'hzsA',
                   'K20933': 'hzsB',
                   'K20932': 'hzsC', }
    for ko in used_ko:
        if ko in ko2gname or ko in rename_dict:
            name = rename_dict.get(ko, '')
            if not name:
                name = ko2gname.get(ko, '')
            koSingle2name[ko] = name
        else:
            r = kegg.get(ko)
            result = [[n for n in _.split(' ') if n] for _ in r.split('\n') if 'NAME' in _]
            name = [_[1] for _ in result][0].strip(',')
            koSingle2name[ko] = name

    with pd.ExcelWriter(join(odir, 'relative_genes_summary.xlsx')) as writer:
        for level in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
            sub_df = locus2info_df.copy()
            sub_df.loc[:, level] = sub_df.loc[:, level].replace('', 'unclassified').fillna('unclassified')
            # total_count = sub_df.loc[:, level].value_counts()
            total_count = sample2info.loc[:, level].replace('', 'unclassified').fillna('unclassified').value_counts()

            collect_dfs = []
            for m in koSingle2name.keys():
                _df = sub_df.loc[sub_df.loc[:, 'ko(single)'] == m, :]

                # count_data = _df.loc[:, level].value_counts()
                count_data, level2snames = count_method(_df, level)
                freq_data = count_data / total_count * 100
                # freq_data = freq_data[freq_data >= 0.6]
                collect_dfs.append(pd.DataFrame(freq_data.values.reshape(-1, 1),
                                                columns=[koSingle2name[m]],
                                                index=freq_data.index))
                summary_df = pd.concat(collect_dfs, axis=1, sort=True)

            summary_df.index = ['%s (%s)' % (_, total_count[_]) for _ in summary_df.index]
            summary_df = summary_df.fillna(0)
            summary_df = summary_df.applymap(lambda x: round(x, 2))
            summary_df = summary_df.reindex(columns=sorted(summary_df.columns))
            summary_df.to_excel(writer, sheet_name=level, index_label=level + ' (in total)')
