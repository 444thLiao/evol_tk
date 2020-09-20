import io
from os.path import expanduser, abspath

import pandas as pd
from Bio import Entrez
from ete3 import NCBITaxa
import xml.etree.ElementTree as ET
from bin.ncbi_convertor import edl
from global_search.thirty_party.metadata_parser import parse_bioproject_xml, parse_biosample_xml

ncbi = NCBITaxa()
taxons = ['superkingdom', 'phylum', 'class',
          'order', 'family', 'genus', 'species']


def tax2tax_info(taxid):
    lineage = ncbi.get_lineage(int(taxid))
    rank = ncbi.get_rank(lineage)
    rank = {v: k for k, v in rank.items()}
    names = ncbi.get_taxid_translator(lineage)
    _dict = {}
    for c in taxons:
        if c in rank:
            _dict[c] = names[rank[c]]
    return _dict


def process_path(path):
    if not '/' in path:
        path = './' + path
    if path.startswith('~'):
        path = expanduser(path)
    if path.startswith('.'):
        path = abspath(path)
    return path


def parse_ipg(t):
    """
    parse result efetched from Identical Protein Groups(ipg)
23582877        INSDC   HQ650005.1      1       452     +       AEG74045.1      ammonia monooxygenase subunit alpha     uncultured bacterium
186872893       INSDC   MG992186.1      1       531     +       AWG96903.1      particulate methane monooxygenase subunit A     uncultured bacterium
    :param t:
    :return:
    """
    bucket = []
    whole_df = pd.read_csv(io.StringIO(t), sep='\t', header=None)
    gb = whole_df.groupby(0)
    all_dfs = [gb.get_group(x) for x in gb.groups]
    for indivi_df in all_dfs:
        indivi_df.index = range(indivi_df.shape[0])
        # aid = indivi_df.iloc[0, 6]
        indivi_df = indivi_df.fillna('')
        pos = [(row[2], row[3], row[4], row[5])
               for row in indivi_df.values]
        gb = [row[10]
              for row in indivi_df.values]
        for row in indivi_df.values:
            bucket.append((row[6], pos, gb))
    return bucket

def parse_xml(xml_data):
    tree = ET.fromstring(xml_data)
    pass

def get_bioproject(bp_list):
    results, failed = edl.esearch(db='bioproject',
                                  ids=bp_list,
                                  result_func=lambda x: Entrez.read(io.StringIO(x))['IdList'])
    all_GI = results[::]
    results, failed = edl.efetch(db='bioproject',
                                 ids=all_GI,
                                 retmode='xml',
                                 retype='xml',
                                 result_func=lambda x: parse_bioproject_xml(x))
    bp2info = {}
    for _ in results:
        if isinstance(_, dict):
            bp2info.update(_)
    return bp2info


def get_biosample(bs_list):
    results, failed = edl.esearch(db='biosample',
                                  ids=bs_list,
                                  result_func=lambda x: Entrez.read(io.StringIO(x))['IdList'])
    all_GI = results[::]
    results, failed = edl.efetch(db='biosample',
                                 ids=all_GI,
                                 retmode='xml',
                                 retype='xml',
                                 result_func=lambda x: parse_biosample_xml(x))
    bs2info = {}
    for _ in results:
        if isinstance(_, dict):
            bs2info.update(_)
    return bs2info
