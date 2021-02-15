import pandas as pd
from bin.ncbi_convertor.func_factory import NCBI_convertor
from bin.ncbi_convertor.toolkit import EntrezDownloader
from bin.ncbi_convertor.pid2bio import genomeID2Bio

tab = pd.read_csv("/home-user/jjtao/Rhizobiales/FLnif-query/gene/query_result/nifH_custom.blast",sep='\t',header=None)
all_nuc_id = tab[0]
all_nuc_id = list(set(all_nuc_id))

all_nuc_id = all_nuc_id[:50] # for testing

edl = EntrezDownloader(
    # An email address. You might get blocked by the NCBI without specifying one.
    # you should use your own account instead of reaching api rate limit.
    # found here: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
    email='l0404th@gmail.com',
    # An API key. You can obtain one by creating an NCBI account. Speeds things up.
    api_key='8ed14220bcca55509d656978cb3f3aa09708',
    num_threads=2,                       # The number of parallel requests to make
    # large threads might raise 'API rate limit exceeded' problems.
    batch_size=50,                        # The number of IDs to fetch per request
    pbar=True                             # Enables a progress bar, requires tqdm package
)
from Bio import Entrez
import io
all_GI = list(unique_id)
results, failed = edl.elink(dbfrom='nuccore',
                                db='assembly',
                            ids=all_GI,
                            result_func=lambda x: Entrez.read(
                                                io.StringIO(x)))
# GI of nuccore to GI of assembly
nid2assembly_dict = {}
for r in results:
    if r['LinkSetDb']:
        nid2assembly_dict[r['IdList'][0]] = r['LinkSetDb'][0]['Link'][0]['Id']

gi_assembly = [_ for _ in nid2assembly_dict.values() if _]
results, failed = edl.esummary(db='assembly',
                                ids=gi_assembly,
                                result_func=lambda x: parse_assembly_xml(x))
# return results is a list of defaultdict.
dbsummary = {}
for _dict in results:
    _dict = dict(_dict)
    for aid, info_dict in _dict.items():
        info_dict['TaxId'] = info_dict['SpeciesTaxid']
        
    dbsummary.update(_dict)


nc = NCBI_convertor(all_nuc_id,'nuccore',given_edl=edl)
nc.get_GI()
GI2nid = {v:k for k,v in nc.GI.items()}
# if an error occured, like "An error occurred, you may see a response text here: {"error":"API rate limit exceeded","
# then some id would be missing 
# you should use nc.get_GI(method='update') to ressuce those missing ids.
all_GI = list(nc.GI.values())
nid2assembly_dict = nc.nid2assembly(all_GI)
# only return those nuccore id with matching assembly id
# And thus, it might empty!!!!!!!!!!!!!!
assembly_GIs = list(nid2assembly_dict.values())
gid2assembly_info, bp2info, bs2info = genomeID2Bio(assembly_GIs)
GI2assembly_id = {_dict["GI"]:aid for aid,_dict in gid2assembly_info.items()}

nid2aid = {GI2nid[nid]:GI2assembly_id.get(aid_GI,'failed') for nid,aid_GI in nid2assembly_dict.items()}




