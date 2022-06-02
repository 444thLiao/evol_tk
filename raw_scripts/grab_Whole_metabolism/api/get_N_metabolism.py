from bioservices.kegg import KEGG
from collections import defaultdict
import pandas as pd
from tqdm import tqdm
import click
import time


def batch_iter(iter, batch_size):
    # generating batch according batch_size
    iter = list(iter)
    n_iter = []
    batch_d = 0
    for batch_u in range(0, len(iter), batch_size):
        if batch_u != 0:
            n_iter.append(iter[batch_d:batch_u])
        batch_d = batch_u
    n_iter.append(iter[batch_d: len(iter) + 1])
    return n_iter


def get_module_info(metabolism_id):
    # N_pathway_id = "ko00910"  # id of N-metabolism

    data = kegg.get(metabolism_id)
    dict_data = kegg.parse(data)

    # get module
    list_module = dict_data.get('MODULE', {})
    module_dict = defaultdict(dict)
    if not metabolism_id.lower().startswith('ko'):
        if not list_module:
            module_ID = 'others'
            module_name = ''
        else:
            module_ID = list(list_module)[0]
            module_name = list(list_module.values())[0]
        # pass a orthology ID
        module_dict[module_ID]['data'] = {}
        module_dict[module_ID]['data']['ORTHOLOGY'] = {metabolism_id: dict_data.get('DEFINITION','')}
        module_dict[module_ID]['data']['NAME'] = [module_name]
    else:
        for m in list_module:
            module_data = kegg.parse(kegg.get(m))
            module_dict[m]['data'] = module_data
            module_dict[m]['metadata'] = list_module[m]
        module_dict['others'] = {}
        module_dict['others']['data'] = {}
        module_dict['others']['data']['ORTHOLOGY'] = dict_data['ORTHOLOGY']
        module_dict['others']['data']['NAME'] = ['Others']
    return module_dict


def assign_ko2info(total_ko2info, ko2info, ko_id):
    if ko_id not in total_ko2info:
        total_ko2info[ko_id]['gene name'] = ';'.join(ko2info.get('NAME', ''))
        total_ko2info[ko_id]['gene definition'] = ko2info.get('DEFINITION', '')
        total_ko2info[ko_id]['reference'] = ';'.join([ref.get('TITLE', '')
                                                      for ref in ko2info.get('REFERENCE', [{}])])
        total_ko2info[ko_id]['num genes'] = len(ko2info.get("GENES", []))
        _cache = list(ko2info.get('MODULE', {}).values())
        if _cache:
            total_ko2info[ko_id]['module'] = ';'.join(_cache)
    else:
        return
    return total_ko2info


def get_locusID(module_dict, locusID_list_output):
    # get orthology/enzyme
    # why don't we use enzyme
    # problematic case:
    # K00372+K00360
    # assimilatory nitrate reductase [EC:1.7.99.-] [RN:R00798]
    ofile = locusID_list_output
    # ofile = '/home-user/thliao/data/metagenomes/N-cycle_locus.list'
    f = open(ofile, 'w')
    locus2info = defaultdict(list)
    total_ko2info = defaultdict(dict)
    for m, mdata in module_dict.items():
        # m: module id ,e.g M00804
        # mdata: dict with "data" & "metadata"
        #   metadata just the name
        #   data contains full informations
        _cache_dict = mdata['data']
        names = [ko_complete.split('+')
                 for ko_complete in _cache_dict['ORTHOLOGY']]
        # names is a list of list
        # each_names is a list of ids. e.g ['K10944','K10945','K10946']
        subunit_ids = [(same_func_unit, subunit_id, full_enzyme_ids)
                       for full_enzyme_ids in names
                       for subunit_id in full_enzyme_ids
                       for same_func_unit in subunit_id.split(',')]
        tqdm.write("processing module %s, it contains %s subunits" % (m, len(subunit_ids)))
        for (same_func_unit, subunit_id, full_enzyme_ids) in tqdm(subunit_ids):
            ko2info = 400
            while isinstance(ko2info, int):
                ko2info = kegg.parse(kegg.get(same_func_unit.strip()))
                time.sleep(0.1)
                # don't make it too fast, or it will block the ip
            new_dict = assign_ko2info(total_ko2info, ko2info, ko_id=same_func_unit)
            if new_dict is not None:
                total_ko2info.update(new_dict)
            org2locus = ko2info['GENES']
            locus_list = [[':'.join([org.lower(), paralog])
                           for paralog in genes.split(' ')]
                          for org, genes in org2locus.items()]
            for locus_paralogs in locus_list:
                for paralog_locus in locus_paralogs:
                    # paralog_locus like e.g. noc:Noc_0892
                    other_paralog_locus = set(locus_paralogs).difference({paralog_locus})
                    other_paralog_locus = ';'.join(list(other_paralog_locus))
                    if '(' in paralog_locus:
                        paralog_locus = paralog_locus.partition('(')[0]
                    if paralog_locus not in locus2info:
                        print(paralog_locus, file=f)
                    locus2info[paralog_locus].append((other_paralog_locus,
                                                      ';'.join(_cache_dict.get('NAME')),
                                                      '+'.join(full_enzyme_ids),
                                                      same_func_unit,
                                                      ))
    f.close()
    return locus2info, total_ko2info


def get_locusDetailedInfo(locus2info):
    genes_df = pd.DataFrame(columns=["locus_name",
                                     "Name",
                                     "source_org",
                                     "definition",
                                     "uniprot ID",
                                     "NCBI ID",
                                     'paralog(within org)',
                                     "module Name",
                                     "Orthology(total)",
                                     "Orthology(single)",
                                     "AA seq", ])
    count_ = 0
    all_locus = list(locus2info.keys())
    pack10_up = batch_iter(all_locus, 10)
    # pack it up, make it 10 times faster
    for bin10 in tqdm(pack10_up):
        bin10 = [_ for _ in bin10 if _]
        info_str = kegg.get('+'.join(bin10))
        for each_str in info_str.split('\nENTRY '):
            if each_str.startswith('ENTRY'):
                info_dict = kegg.parse(each_str)
            else:
                info_dict = kegg.parse('ENTRY ' + each_str)
            locus = info_dict.get('ENTRY', 'unknown').split(' ')[0]
            try:
                locus = [_ for _ in bin10 if locus.lower() in _.lower()][0]
            except:
                locus = info_dict.get('NAME', [])[0]
                locus = [_ for _ in bin10 if locus.lower() in _.lower()][0]
                import pdb;pdb.set_trace()
            if locus in genes_df.iloc[:, 0]:
                continue

            args = locus2info[locus][0]
            other_paralog_locus, module_name, orthology_total, orthology_single = args
            gene_name = ';'.join(info_dict.get('NAME', ['']))
            definition = info_dict['DEFINITION']
            source_organism = info_dict["ORGANISM"]
            NCBI_refID = info_dict.get("DBLINKS", {}).get("NCBI-ProteinID", None)
            uniprot_refID = info_dict.get("DBLINKS", {}).get("UniProt", None)
            AA_seq = info_dict["AASEQ"].replace(' ', '')
            genes_df.loc[count_, :] = [locus,
                                       gene_name,
                                       source_organism,
                                       definition,
                                       uniprot_refID,
                                       NCBI_refID,
                                       other_paralog_locus,
                                       module_name,
                                       orthology_total,
                                       orthology_single,
                                       AA_seq]
            count_ += 1
    return genes_df


def main(metabolism_id, locus_id_list, extra_ko):
    tqdm.write("getting module information")
    module_dict = get_module_info(metabolism_id)
    if extra_ko is not None:
        module_dict['extra'] = {}
        module_dict['extra']['data'] = {}
        module_dict['extra']['data']['ORTHOLOGY'] = extra_ko
        module_dict['extra']['data']['NAME'] = ['Others']
    tqdm.write("For each module, getting it locus/gene ID and sequence")
    locus2info, all_ko2info = get_locusID(module_dict, locus_id_list)
    ko_df = pd.DataFrame.from_dict(all_ko2info, orient='index')
    return ko_df, locus2info
    # tqdm.write("get locus infomation from database, and also get sequence")
    # genes_df = get_locusDetailedInfo(locus2info)
    # return genes_df, ko_df


@click.command()
@click.option("-ko_id", "koID")
@click.option("-locusID", "locusID_out", help="output locus ID list")
@click.option("-koDF", "koDF_out", help="output table of relative ko number")
@click.option("-locusDF", "locusDF_out", help="output table of relative locus ID and its sequence")
@click.option("-drop_ko", "removed_ko", help="file indicated ko doesn't need", default=None)
@click.option("-extra_ko", "extra_ko", help="file indicated extra ko", default=None)
def cli(koID, locusID_out, koDF_out, locusDF_out, removed_ko, extra_ko):
    if extra_ko is not None:
        extra_ko = set([_ for _ in open(extra_ko).read().split('\n') if _])
    ko_df, locus2info = main(koID, locusID_out, extra_ko)
    ko_df.to_csv(koDF_out, sep='\t', index=1, index_label="KO number")

    if removed_ko is not None:
        removed_ko = set([_ for _ in open(removed_ko).read().split('\n') if _])
        locus2info = {k: v
                      for k, v in locus2info.items()
                      if v[-1][-1] not in removed_ko}
    genes_df = get_locusDetailedInfo(locus2info)
    genes_df.loc[:, 'KO name'] = [ko_df.loc[ko, 'gene name']
                                  for ko in genes_df.loc[:, 'Orthology(single)']]
    genes_df.to_csv(locusDF_out, sep='\t', index=0)

if __name__ == '__main__':
    kegg = KEGG()
    # init a engine to fetch kegg id
    cli()
    # python3 get_N_metabolism.py -locusID ./Nitrogen_cycle_locus.list -koDF ./ko_info.csv -locusDF ./locus_info.csv -ko_id ko00910
    # metabolism_id = "ko00910"  # id of N-metabolism
    # locus_id_list = '/home-user/thliao/data/metagenomes/N-cycle_locus.list'
    # genes_df, ko_df = main(metabolism_id, locus_id_list)
    # ko_df.to_csv("/home-user/thliao/data/metagenomes/KO_info.tsv", sep='\t', index=1,index_label="KO number")
    # genes_df.to_csv("/home-user/thliao/data/metagenomes/N-relative_genes.tsv", sep='\t', index=0)
    # locus_tag could not used as index, because it has duplicated
