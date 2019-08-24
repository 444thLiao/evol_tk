from bioservices.kegg import KEGG
from collections import defaultdict
import pandas as pd
from tqdm import tqdm

k = KEGG()


def get_module_info(metabolism_id):
    # N_pathway_id = "ko00910"  # id of N-metabolism

    data = k.get(metabolism_id)
    dict_data = k.parse(data)

    # get module
    list_module = dict_data['MODULE']
    module_dict = defaultdict(dict)
    for m in list_module:
        module_data = k.parse(k.get(m))
        module_dict[m]['data'] = module_data
        module_dict[m]['metadata'] = list_module[m]
    return module_dict


def get_locusID(module_dict, locusID_list_output):
    # get orthology/enzyme
    # why don't we use enzyme
    # problematic case:
    # K00372+K00360
    # assimilatory nitrate reductase [EC:1.7.99.-] [RN:R00798]

    ofile = locusID_list_output
    # ofile = '/home-user/thliao/data/metagenomes/N-cycle_locus.list'
    f = open(ofile, 'w')
    m2orthology = defaultdict(dict)
    locus2info = defaultdict(list)
    for m, mdata in module_dict.items():
        # module id
        _cache_dict = mdata['data']
        names = [_.split('+') for _ in _cache_dict['ORTHOLOGY']]
        m2orthology[m]['metadata'] = _cache_dict['ORTHOLOGY']
        for each_names in tqdm(names):
            m2orthology[m]['orthology'] = {}
            for name in each_names:
                # othology id/ gene id but separate with '+' which means subunit
                # maybe each_names contains only one, because no subunit
                for same_func_name in name.split(','):
                    org2locus = k.parse(k.get(same_func_name.strip()))['GENES']
                    locus_list = [[':'.join([k.lower(), paralog]) for paralog in v.split(' ')]
                                  for k, v in org2locus.items()]
                    for locus_paralogs in locus_list:

                        for paralog_locus in locus_paralogs:
                            other_paralog_locus = set(locus_paralogs).difference({paralog_locus})
                            other_paralog_locus = ';'.join(list(other_paralog_locus))
                            if '(' in paralog_locus:
                                paralog_locus = paralog_locus.partition('(')[0]
                            if paralog_locus not in locus2info:
                                print(paralog_locus, file=f)
                            # if paralog_locus in locus2info:
                            #     print('!!!', locus2info[paralog_locus], [other_paralog_locus,
                            #                                              ';'.join(_cache_dict.get('NAME')),
                            #                                              '+'.join(each_names),
                            #                                              same_func_name,
                            #                                              ])
                            locus2info[paralog_locus].append((other_paralog_locus,
                                                              ';'.join(_cache_dict.get('NAME')),
                                                              '+'.join(each_names),
                                                              same_func_name,
                                                              ))
    f.close()
    return locus2info


def get_locusDetailedInfo(locus2info):
    genes_df = pd.DataFrame(columns=["Name",
                                     "source_org",
                                     "definition",
                                     "uniprot ID",
                                     "NCBI ID",
                                     'paralog(within org)',
                                     "module Name",
                                     "Orthology(total)",
                                     "Orthology(single)",
                                     "AA seq", ])
    for locus, locus_info_list in tqdm(locus2info.items(), total=len(locus2info)):
        if locus in genes_df.index:
            continue
        info = k.parse(k.get(locus))
        for (other_paralog_locus,
             module_name,
             orthology_total,
             orthology_single) in locus_info_list:

            if type(info) == dict:
                gene_name = ';'.join(info.get('NAME', ['']))
                definition = info['DEFINITION']
                source_organism = info["ORGANISM"]
                NCBI_refID = info.get("DBLINKS", {}).get("NCBI-ProteinID", None)
                uniprot_refID = info.get("DBLINKS", {}).get("UniProt", None)
                AA_seq = info["AASEQ"].replace(' ', '')
                genes_df.loc[locus, :] = [gene_name,
                                          source_organism,
                                          definition,
                                          uniprot_refID,
                                          NCBI_refID,
                                          other_paralog_locus,
                                          module_name,
                                          orthology_total,
                                          orthology_single,
                                          AA_seq]
    return genes_df


def main(metabolism_id):
    module_dict = get_module_info(metabolism_id)
    locus2info = get_locusID(module_dict)
    genes_df = get_locusDetailedInfo(locus2info)
    return genes_df


if __name__ == '__main__':
    metabolism_id = "ko00910"  # id of N-metabolism
    genes_df = main(metabolism_id)
