import multiprocessing as mp
import os
from collections import Counter
from os.path import dirname, abspath, exists, join
from subprocess import check_call
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
from bioservices import KEGG


def output_seq_from_df(df, oseq):
    with open(oseq, 'w') as f1:
        for _, row in df.iterrows():
            aa_seq = row['AA seq']
            name = row["locus_name"]
            f1.write(f'>{name}\n')
            f1.write(f"{aa_seq}\n")


def run_cmd(cmd):
    check_call(cmd, shell=True)


used_ko = []
removed_ko = [_ for _ in open('./removed_ko.txt').read().split('\n')
              if _]
used_ko += [_.split('\t')[0] for _ in open('./ko_info.csv').read().split('\n')
            if _ and _.split('\t')[0] != 'KO number']
used_ko = [_ for _ in used_ko if _ not in removed_ko]

# relative to /home-user/thliao/data/metagenomes/manually_output
# input file
locus2info = './locus_info.csv'  # genes from kegg database
sample2locus = '../concat_all/sample2infos.tsv'  # collected metadata
manually_info = '../manually_curated_N_cycle_genes.xlsx'  # manually curated genes with ko info? or not

ko2info_file = '../KO_info.xlsx'
target_fa = '../concat_all/all_protein.faa'

# output file
odir = '.'
oseq = join(odir, 'first_extract_seq.faa')
o1_tab = join(odir, 'first_diamond.out')
intermedia_faa = join(odir, 'first_extract_seq.faa')
o2_tab = join(odir, 'first_extract_seq_KOfam.out')
final_faa = join(odir, 'confirmed.faa')
final_tsv = join(odir, 'confirmed_locus2info.tsv')


def main(locus2info, sample2locus, target_fa, oseq):
    # step1. extract genes info
    df = pd.read_csv(locus2info, sep='\t')
    odir = dirname(oseq)
    if not exists(odir):
        os.makedirs(odir)
    output_seq_from_df(df, oseq)
    dbname = abspath(oseq).rpartition('.')[0]

    # build index with step1 output
    run_cmd(f'diamond makedb --in {oseq} --db {dbname}')
    # perform first blastp
    run_cmd(f"diamond blastp -q {target_fa} -o {o1_tab} -d {dbname} -p 0 -b 5 -c 2")

    # reannotate from gene info file, especially the KO and its name
    pre_df = pd.read_csv(f'{o1_tab}', sep='\t', header=None)
    subject_info_df = pd.read_excel(manually_info)
    subject_info_df = subject_info_df.set_index('AA accession')
    order_df = subject_info_df.reindex(pre_df.loc[:, 1])
    pre_df.loc[:, 'cover ratio'] = pre_df.loc[:, 3].values / order_df.loc[:, 'AA sequence(seq)'].str.len().values
    pre_df.loc[:, 'KO'] = order_df.loc[:, 'ko'].values
    pre_df.loc[:, 'KO name'] = order_df.loc[:, 'gene name'].values
    # subject_info_df = subject_info_df.set_index('locus name')
    # query_df = query_df.drop_duplicates('locus_name')
    # query_df = query_df.set_index('locus_name')
    # get_df = query_df.reindex(pre_df.loc[:, 1])
    # pre_df.loc[:, 'cover ratio'] = pre_df.loc[:, 3].values / get_df.loc[:, 'AA seq'].str.len().values
    # pre_df.loc[:, 'KO'] = get_df.loc[:, 'Orthology(single)'].values
    # pre_df.loc[:, 'KO name'] = get_df.loc[:, 'KO name'].values

    # considerate before and after locus, which being cutted
    _pre_df = pre_df.drop_duplicates(0)
    all_locus_ed = set(_pre_df.loc[:, 0].unique())
    counted_locus = []
    for _, row in tqdm(_pre_df.iterrows(), total=_pre_df.shape[0]):
        locus = row[0]
        n = row['KO name']
        num = int(locus[-2:])
        b_locus, after_locus = locus[:-2] + str(int(num - 1)), locus[:-2] + str(int(num + 1))
        if b_locus in all_locus_ed:
            _sub_df = _pre_df.loc[_pre_df.loc[:, 0] == b_locus, 'KO name'].values[0]
            if _sub_df == n:
                counted_locus.append(tuple(sorted([locus, b_locus]) + [n]))
                # print(locus,b_locus,n)
        if after_locus in all_locus_ed:
            _sub_df = _pre_df.loc[_pre_df.loc[:, 0] == after_locus, 'KO name'].values[0]
            if _sub_df == n:
                counted_locus.append(tuple(sorted([locus, after_locus]) + [n]))
                # print(locus,after_locus,n)

    tmp_df = pd.read_csv(f'{o1_tab}', sep='\t', header=None)
    records = SeqIO.parse(f'{target_fa}', format='fasta')
    used_gids = set(tmp_df.iloc[:, 0])
    collcect_records = []
    for record in tqdm(records):
        if record.id in used_gids:
            collcect_records.append(record)
    with open(f'{intermedia_faa}', 'w') as f1:
        SeqIO.write(collcect_records, f1, format='fasta-2line')

    # step2: using more and well-annotated kegg database to annotated, for removing false-positive
    run_cmd(f"/home-user/thliao/software/kofamscan/exec_annotation {intermedia_faa} -o {o2_tab} --cpu 2 -f mapper-one-line")
    # extract step2 output info into a dictionary
    kofamscan_str = open(f'{o2_tab}', 'r').read().split('\n')
    ko2result = dict([(row.split('\t')[0],
                       row.split('\t')[1:])
                      for row in kofamscan_str
                      ])

    # step3: init a bucket to collect all required info
    confirmed_seq = defaultdict(dict)
    # prepare annotation
    ko_info = pd.read_excel(ko2info_file)
    ko2name = dict(zip(ko_info.iloc[:, 0].values, ko_info.iloc[:, 1].values))
    ko2name.update(dict([('K20934', 'hzsA'),
                         ('K20933', 'hzsB'),
                         ('K20932', 'hzsC'),
                         ('K12265', 'norW'),
                         ('K12264', 'norV'),
                         ('K05916', 'hmp, YHB1')]))

    # pre_ko2name = dict(zip(subject_info_df.loc[:, 'ko'].values,
    #                        subject_info_df.loc[:, 'gene name'].values))
    # pre_ko2name.pop('nan')  # nan maybe multiple, could not overwrite by one
    # ko2name.update(pre_ko2name)
    def annotate_KO(locus):
        collect_seq = defaultdict(dict)
        sub_df = pre_df.loc[pre_df.loc[:, 0] == locus, :]
        _df = sub_df.loc[(sub_df.loc[:, 'cover ratio'] >= 0.6) & (sub_df.loc[:, 10] <= 1e-7), :]
        ko1_list = []
        pre_ko2name = {}
        ko_output = []
        existing_names = []
        if _df.shape[0] != 0:
            existing_names = list(_df['KO name'].fillna('missing'))
            ko1_list = list(_df["KO"].fillna('missing'))
            ko_output = ['%s:%s' % (k, n)
                         for k, n in zip(ko1_list, existing_names)]
            pre_ko2name = dict(zip(ko1_list, existing_names))
            if 'missing' in pre_ko2name:
                pre_ko2name.pop('missing')  # nan maybe multiple, could not overwrite by one
        annotated_KOs_list = ko2result.get(locus, [])
        annotated_KOs_set = set(annotated_KOs_list)
        if not annotated_KOs_list and not ko1_list:
            return collect_seq
        elif annotated_KOs_list:
            top1_annotated = annotated_KOs_list[0]
            if top1_annotated in set(ko1_list):
                # in pre-annotated set, maybe totally same
                collect_seq[locus]['KO'] = top1_annotated
                collect_seq[locus]['KO name'] = pre_ko2name.get(top1_annotated, '')
                if len(set(ko1_list)) == 1:
                    collect_seq[locus]['within paralog'] = 'NO'
                else:
                    collect_seq[locus]['within paralog'] = 'YES'
                collect_seq[locus]['outside paralog'] = 'NO'
            elif top1_annotated in ko2name:
                # not in pre-annotated set, but in total set
                collect_seq[locus]['KO'] = top1_annotated
                collect_seq[locus]['KO name'] = ko2name[top1_annotated]
                collect_seq[locus]['within paralog'] = 'YES'
                collect_seq[locus]['outside paralog'] = 'NO'
            else:
                # not in .., and not in .., total new or missing ko.
                collect_seq[locus]['KO'] = top1_annotated
                collect_seq[locus]['KO name'] = "new one"
                collect_seq[locus]['within paralog'] = 'NO'
                collect_seq[locus]['outside paralog'] = 'YES'
            _ko_output = [o
                          for o in ko_output
                          if o.split(':')[0] not in annotated_KOs_set]
            if _ko_output:
                collect_seq[locus]['likely KO'] = [_ko_output[0]]

        elif ko1_list:
            top1_pre = ko1_list[0]
            top1_name = existing_names[0]
            collect_seq[locus]['KO'] = top1_pre
            collect_seq[locus]['KO name'] = top1_name
            collect_seq[locus]['within paralog'] = 'NO'
            collect_seq[locus]['outside paralog'] = 'NO'
        return collect_seq

    unique_locus_name_indf = list(set(pre_df.loc[:, 0]))
    with mp.Pool(processes=64) as tp:
        for r in tqdm(tp.imap(annotate_KO, unique_locus_name_indf),
                      total=len(unique_locus_name_indf)):
            confirmed_seq.update(r)
    # collect all confirmed id
    collect_IDs = [_
                   for _, v in confirmed_seq.items()
                   if (v.get('likely KO')) or
                   v.get('outside paralog') != 'YES']

    collect_IDs_set = set(collect_IDs)
    print("contains %s confirmed sequence" % len(collect_IDs_set))
    records = SeqIO.parse(f'{intermedia_faa}', format='fasta')
    collect_reads = [_
                     for _ in records
                     if _.id in collect_IDs_set]

    with open(f'{final_faa}', 'w') as f1:
        SeqIO.write(collect_reads, f1, format='fasta-2line')

    sample2info = pd.read_csv(sample2locus, sep='\t', header=0, index_col=1)
    locus2info_df = pd.DataFrame(columns=["locus_prefix",
                                          'sample name',
                                          'source project',
                                          'ko(single)',
                                          'Gene name(N metabolism)',
                                          'within paralog',
                                          'paralog KO',
                                          'paralog KO name',
                                          ])
    new_dict = defaultdict(dict)
    for locus in tqdm(collect_IDs_set):
        _dict = confirmed_seq[locus]
        ko = _dict['KO']
        ko_name = _dict['KO name']
        locus_prefix = locus.split('_')[0]
        sname, sproject = sample2info.loc[locus_prefix, ['sample_name', 'source']].values
        wp, wpk = _dict.get('within paralog', 'NO'), _dict.get('likely KO', '')
        for _, v in zip(locus2info_df.columns, [locus_prefix, sname, sproject, ko, ko_name, wp,
                                                ';'.join([_.split(':')[0] for _ in wpk]),
                                                ';'.join([_.split(':')[1] for _ in wpk])]):
            new_dict[locus][_] = v
    locus2info_df = pd.DataFrame.from_dict(new_dict, orient='index')
    _sub_df = locus2info_df.loc[locus2info_df.loc[:, 'Gene name(N metabolism)'] == 'new one', ['paralog KO',
                                                                                               'paralog KO name']]
    locus2info_df.loc[locus2info_df.loc[:, 'Gene name(N metabolism)'] == 'new one',
                      ['ko(single)',
                       'Gene name(N metabolism)']] = _sub_df.values
    locus2info_df = locus2info_df.loc[locus2info_df.loc[:, 'ko(single)'] != 'K00373', :]
    # drop k00373

    locus2info_df.loc[:, 'Gene name(N metabolism)'] = list(map(lambda x: x[1] if x[0] == 'missing' else ko2name[x[0]],
                                                               locus2info_df.loc[:, ['ko(single)',
                                                                                     'Gene name(N metabolism)']].values))
    locus2info_df.loc[locus2info_df.loc[:, 'paralog KO name'] == 'norZ',
                      'Gene name(N metabolism)'] = 'norZ'
    # detect duplicated one, remove smallest one?
    _count = []
    for _ in set(counted_locus):
        if _[0] in locus2info_df.index and _[1] in locus2info_df.index:
            if _[2] in locus2info_df.loc[_[0], ['Gene name(N metabolism)', 'paralog KO name']].values:
                _count.append(_)
    for d_tuple in tqdm(_count):
        l1, l2 = d_tuple[:2]
        dropped_locus = sorted([l1, l2], key=lambda x: pre_df.loc[pre_df.loc[:, 0] == x, 3].values[0])[0]
        locus2info_df = locus2info_df.drop(dropped_locus)
    t = ['CIMKLAMG_01008', 'BLHCAHHO_01887', 'KAPHPNMJ_01950',
         'BLHCAHHO_02203', 'KAPHPNMJ_00555', 'GFGMIJHE_02191',
         'DDDAEOCG_00928', 'FPKBIPEM_01264', 'OKOMEFJD_00369',
         'DFJGGOCD_00279', 'FNMJIPCB_00158', 'KKEGIMMP_00878',
         'GLELHIOC_00437', 'CNAAHEGP_00473', 'ILGCLMOP_00797',
         'OMKEHILC_01648', 'IEGHJKID_00458', 'GFGMIJHE_01595',
         'DDDAEOCG_01667', 'KKHHJNPO_01081', 'LAIKLKFG_01802',
         'KDLDFFLJ_00880', 'AMMCEENH_00385', 'OKOMEFJD_00123',
         'EHGCCCFP_01499', 'FDGGIJME_00679', 'FPKBIPEM_01592',
         'KFHEOECM_00769', 'GLELHIOC_01843', 'ILGCLMOP_00106',
         'OMKEHILC_01002', 'IEGHJKID_01654', 'CNAAHEGP_00879',
         'KKEGIMMP_00079', 'KFHEOECM_01249', 'BLHCAHHO_00728',
         'KAPHPNMJ_01068', 'LAIKLKFG_00187', 'FPKBIPEM_01120',
         'ILGCLMOP_00394', 'OMKEHILC_01706', 'IEGHJKID_01360',
         'GLELHIOC_00795', 'CNAAHEGP_00245', 'BLHCAHHO_02199',
         'KAPHPNMJ_00009', 'GFGMIJHE_02168', 'DDDAEOCG_00821',
         'KFHEOECM_01063', 'ILGCLMOP_00801', 'OMKEHILC_01644',
         'IEGHJKID_00454', 'GLELHIOC_00441', 'CNAAHEGP_00469',
         'FPKBIPEM_01268', 'ENIBLNEE_00528', 'FNMJIPCB_01245',
         'EJGMDGIB_00419', 'DFJGGOCD_00276', 'FDGGIJME_00006',
         'OKOMEFJD_00372'],
    locus2info_df.loc[t, 'Gene name(N metabolism)'] = 'hdh'
    locus2info_df.loc[t, 'ko(single)'] = 'K20935'
    #
    # for _, row in locus2info_df.iterrows():
    #     ko = row['ko(single)']
    #     name = row['Gene name(N metabolism)']
    # if name == 'unknown' and ko in ko2name:
    #     locus2info_df.at[_, 'Gene name(N metabolism)'] = ko2name[ko]

    order_sample2info = sample2info.copy()
    order_sample2info = order_sample2info.reindex([_
                                                   for _ in locus2info_df.locus_prefix])
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
         'species(from metadata)']
    locus2info_df = locus2info_df.reindex(columns=list(locus2info_df.columns) + c)
    locus2info_df.loc[:, c] = order_sample2info.loc[:, c].values
    locus2info_df.to_csv(final_tsv, sep='\t', index=1)
    ############################################################
    # concat manually df
    previous_result = pd.read_csv('../new_grab/output/confirmed_locus2info.tsv', sep='\t', index_col=0)
    previous_result.columns = locus2info_df.columns
    _sub_df = previous_result.loc[previous_result.index.difference(locus2info_df.index), :]
    final_df = pd.concat([locus2info_df, _sub_df], axis=0)
    final_df.to_csv('./concated_all.csv', sep='\t', index=1)

    locus2info_df = final_df.copy()
    order_columns = ['anfG', 'nifD', 'nifK', 'nifH', 'vnfG', 'vnfD', 'vnfK', 'vnfH',
                     'pmoA-amoA', 'pmoB-amoB', 'pmoC-amoC', 'hao', 'nirK',
                     'narG, narZ, nxrA', 'narH, narY, nxrB', 'pmoA-amoA', 'pmoB-amoB',
                     'pmoC-amoC', 'hao', 'hmp, YHB1', 'nasA', 'nasB',
                     'narG, narZ, nxrA', 'narH, narY, nxrB', 'narI, narV', 'napA',
                     'napB', 'nirK', 'nirS', 'norV', 'norW', 'hcp', 'norB', 'norC', 'norZ',
                     'nosZ', 'narG, narZ, nxrA', 'narH, narY, nxrB', 'narI, narV',
                     'napA', 'napB', 'nrfA', 'nrfH', 'εHao (HaoA)', 'haoC', 'ONR', 'OTR', 'nirK', 'nirS',
                     'hzsA', 'hzsB', 'hzsC', 'hdh', 'hao', 'nasA', 'nasB', 'narB',
                     'nirA', 'nirB', 'nirD']
    with pd.ExcelWriter(join(odir, 'relative_genes_summary(tax).xlsx')) as writer:
        for level in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
            sub_df = locus2info_df.copy()
            sub_df.loc[:, level] = sub_df.loc[:, level].replace('', 'unclassified').fillna('unclassified')
            # total_count = sub_df.loc[:, level].value_counts()
            total_count = sample2info.loc[:, level].replace('', 'unclassified').fillna('unclassified').value_counts()

            collect_dfs = []
            for genename in sub_df.loc[:, 'Gene name(N metabolism)'].unique():
                _df = sub_df.loc[sub_df.loc[:, 'Gene name(N metabolism)'] == genename, :]
                _df = _df.drop_duplicates('locus_prefix')
                count_data = _df.loc[:, level].value_counts()
                # count_data, level2snames = count_method(_df, level)
                freq_data = count_data / total_count * 100
                # freq_data = freq_data[freq_data >= 0.6]
                collect_dfs.append(pd.DataFrame(freq_data.values.reshape(-1, 1),
                                                columns=[genename],
                                                index=freq_data.index))
                summary_df = pd.concat(collect_dfs, axis=1, sort=True)

            summary_df.index = ['%s (%s)' % (_, total_count[_]) for _ in summary_df.index]
            summary_df = summary_df.fillna(0)
            summary_df = summary_df.applymap(lambda x: round(x, 2))
            summary_df = summary_df.reindex(columns=order_columns)
            summary_df = summary_df.reindex(sorted(summary_df.index,
                                                   key=lambda x: int(x.split('(')[1].strip(')')),
                                                   reverse=True))
            summary_df = summary_df.T
            summary_df.to_excel(writer, sheet_name=level, index_label=level + ' (in total)')

    project2env = {'19_Stewart': 'rumen',
                   '17_lee': 'Fecal',
                   '15_brown': 'subsurface aquifer',
                   '16_Anantharaman': 'subsurface aquifer',
                   '16_haroon': 'marine',
                   '17_jungbluth': 'hydrothermal fluid',
                   '17_Tully': 'marine',
                   '17_parks': 'public data',
                   '18_Tully': 'subsurface aquifer',
                   '18_probst': 'subsurface aquifer',
                   '17_Hernsdorf': 'subsurface aquifer',
                   '18_Woodcroft': 'permafrost',
                   '18_crits': 'soil',
                   '18_Vavourakis': 'hypersaline',
                   '19_Rinke': 'marine',
                   '19_pedron': 'spring',
                   '18_Stewart': 'rumen',
                   '18_Delmont': 'marine'}
    env_list = set(project2env.values())

    with pd.ExcelWriter(join(odir, 'relative_genes_summary(env).xlsx')) as writer:
        total_collect = []
        for env in env_list:
            projects = [k for k, v in project2env.items() if v == env]
            sub_df = locus2info_df.copy()
            sub_df = sub_df.loc[sub_df.loc[:, 'source project'].isin(projects), :]
            total_count = sample2info.loc[sample2info.source.isin(projects), :].shape[0]

            collect_dfs = []
            for m in ko2name.keys():
                _df = sub_df.loc[sub_df.loc[:, 'ko(single)'] == m, :]
                _df = _df.drop_duplicates('locus_prefix')
                count_data = _df.loc[:, level].shape[0]
                # count_data, level2snames = count_method(_df, level)
                freq_data = count_data / total_count * 100
                # freq_data = freq_data[freq_data >= 0.6]
                collect_dfs.append(pd.DataFrame(freq_data,
                                                columns=[koSingle2name[m]],
                                                index=['%s (%s)' % (env, total_count)]))
                summary_df = pd.concat(collect_dfs, axis=1, sort=True)
            total_collect.append(summary_df)
        summary_df = pd.concat(total_collect)
        summary_df = summary_df.fillna(0)
        summary_df = summary_df.applymap(lambda x: round(x, 2))
        summary_df = summary_df.reindex(columns=order_columns)
        summary_df.to_excel(writer, index_label='ENV (in total)')

    ############################################################
    DB_locus2info = pd.read_csv('./locus_info.csv', sep='\t', index_col=0)
    DB_locus2info = DB_locus2info.drop_duplicates(['source_org', 'KO name'])
    source_org_list = [_.split(' ')[0] for _ in DB_locus2info.loc[:, 'source_org']]
    from io import StringIO
    db_ref = pd.read_csv(StringIO(kegg.list('organism')), sep='\t', header=None)
    db_ref = db_ref.set_index(1)
    db_ref = db_ref.reindex(source_org_list)
    DB_locus2info.loc[:, 'lineage'] = db_ref.loc[:, 3].values

    gname2lineage = dict()
    for gname in DB_locus2info.loc[:, 'KO name'].unique():
        sub_df = DB_locus2info.loc[DB_locus2info.loc[:, 'KO name'] == gname, 'lineage']
        gname2lineage[gname] = list(set(sub_df))
