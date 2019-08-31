from bioservices.kegg import KEGG
from collections import defaultdict
import pandas as pd
from tqdm import tqdm
import click

def get_module_info(metabolism_id):
    # N_pathway_id = "ko00910"  # id of N-metabolism

    data = kegg.get(metabolism_id)
    dict_data = kegg.parse(data)

    # get module
    list_module = dict_data['MODULE']
    module_dict = defaultdict(dict)
    for m in list_module:
        module_data = kegg.parse(kegg.get(m))
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
    locus2info = defaultdict(list)
    for m, mdata in module_dict.items():
        # m: module id ,e.g M00804
        # mdata: dict with "data" & "metadata"
        #   metadata just the name
        #   data contains full informations
        _cache_dict = mdata['data']
        names = [_.split('+')
                 for _ in _cache_dict['ORTHOLOGY']]
        # names is a list of list
        # each_names is a list of ids. e.g ['K10944','K10945','K10946']
        subunit_ids = [(same_func_unit, subunit_id, full_enzyme_ids)
                       for full_enzyme_ids in names
                       for subunit_id in full_enzyme_ids
                       for same_func_unit in subunit_id.split(',')]

        for (same_func_unit, subunit_id, full_enzyme_ids) in tqdm(subunit_ids):
            org2locus = kegg.parse(kegg.get(same_func_unit.strip()))['GENES']
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
    return locus2info


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
    for locus, locus_info_list in tqdm(locus2info.items(),
                                       total=len(locus2info)):
        if locus in genes_df.index:
            continue
        info_dict = kegg.parse(kegg.get(locus))
        for (other_paralog_locus,
             module_name,
             orthology_total,
             orthology_single) in locus_info_list:

            if type(info_dict) == dict:
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



def main(metabolism_id, locus_id_list):
    module_dict = get_module_info(metabolism_id)
    locus2info = get_locusID(module_dict, locus_id_list)
    genes_df = get_locusDetailedInfo(locus2info)
    return genes_df


@click.command()
@click.argument("ko_id",help='id of a metabolism. e.g for nitrogen metabolism is ko00910 ')
@click.option('-o','ofile')
@click.option("-o_locus","olocus_file")
def cli(ko_id):
    main

if __name__ == '__main__':
    kegg = KEGG()
    metabolism_id = "ko00910"  # id of N-metabolism
    locus_id_list = '/home-user/thliao/data/metagenomes/N-cycle_locus.list'
    genes_df = main(metabolism_id, locus_id_list)
    genes_df.to_csv("/home-user/thliao/data/metagenomes/N-relative_genes.tsv", sep='\t', index=0)
    # locus_tag could not used as index, because it has duplicated
