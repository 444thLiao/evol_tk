import io,os
from collections import defaultdict

from tqdm import tqdm
from Bio import SeqIO, Entrez
from bin.ncbi_convertor.toolkit import (
    edl,
    access_intermedia,
    tax2tax_info,
    parse_ipg,
    get_GI,
    unpack_gb
)
from api_tools.third_party import parse_assembly_xml
import Bio


tax_convertable_dbs = ["protein", "assembly", "nuccore",'nucleotide']
batch_return_dbs = []
single_return_dbs = ["protein", "nuccore",'nucleotide']
db_without_GI  = ["protein",'nuccore']



def eread(x):
    if float(Bio.__version__)>=1.77:
        return Entrez.read(io.BytesIO(x.encode()))
    else:
        return Entrez.read(io.StringIO(x))

def read_efetch(t):
    tmp = {}
    records = list(SeqIO.parse(io.StringIO(t),'genbank'))
    for r in records:
        tmp[r.id] = unpack_gb(r)
    return [tmp]
    
class NCBI_convertor:
    "prototype of function which convert ID to GI"

    def __init__(self, IDs, db, given_edl=edl):
        """
        :param IDs: The requested ID
        :param db: Which database of those ID belong to
        :return:
        """
        self.origin_ids = IDs
        self.dbname = db
        self.GI = None
        self.tids = {}
        self.dbsummary = {}
        self.edl = given_edl
        self.validate_methods = ['get','update']
        if self.edl.email is None:
            print("EDL you pass without an API key. If you pass a key, it could speed up 3times at most. You could use API key by set the environment parameters 'EMAIL' and 'EKEY', and then run it in new terminal. Or use the function refresh_key to update the email and api_key. " )
    
    def refresh_key(self,email=None,api_key=None):
        if email is None:
            email = os.getenv('EMAIL')
            api_key = os.getenv("EKEY")
        self.edl.email = email
        self.edl.api_key = api_key
        
    def load_from_file(self, infile):
        # todo:
        pass

    def print_available_dbs(self):
        print("")

    def get_biosample(self):
        pass

    def get_bioproject(self):
        pass
    def check_methods(self,m):
        if m not in self.validate_methods:
            raise IOError(f'method {m} is not in the validated_methods')
        
    def check_cache(self, suffix, redo):
        _cache = access_intermedia(self.origin_ids, suffix=suffix, redo=redo)
        if _cache is not None:
            id2other_IDs = _cache
            if "2GI" in suffix:
                self.GI = id2other_IDs
            elif "" in suffix:
                pass
        else:
            print("No cache found.")

    def get_seq(self, batch_size=5, method="get",preset='seq'):
        # for protein
        self.check_methods(method)
        
        ids = self.origin_ids
        if self.dbname not in ["nucleotide", "protein", 'nuccore']:
            raise SyntaxError
        if preset == 'seq':
            retype='fasta'
            format='fasta'
        elif preset=='genbank':
            retype='gb'
            format='genbank'
        results, failed = self.edl.efetch(
            db=self.dbname,
            ids=ids,
            result_func=lambda x: list(SeqIO.parse(io.StringIO(x), format)),
            batch_size=batch_size,
            retype=retype,
            retmode="text",
        )
        return results

    def get_GI(self, num_retry=5, method="get"):
        self.check_methods(method)
        if self.dbname in db_without_GI:
            return {k:k for k in self.origin_ids}
        # if it is the db that needed to be retrieved GIs
        if (self.GI is not None) and (method != "update"):
            return self.GI
        if method == "update" and self.GI is not None:
            ids = [k for k, v in self.GI.items() if not v]
        elif method == "get" or self.GI is not None:
            print(f"original values have {len(self.origin_ids)} but returns pre-existed {len(self.GI)}")
            return self.GI
        elif method == "get" or self.GI is None:
            ids = self.origin_ids
        else:
            raise IOError('unexpected if/else')
        
        id2gi, failed = get_GI(ids, self.dbname, self.edl)
        # stodge the result into intermedia file for second access.
        process_name = f"{self.dbname}2GI"
        access_intermedia(id2gi, suffix=process_name)

        if self.GI is None:
            self.GI = id2gi
        else:
            self.GI.update(id2gi)
        self.failed_ids = failed
        return self.GI

    def get_db_summary(self, all_GI=None, method="update"):
        self.check_methods(method)
        if all_GI is None:
            all_GI = list(self.get_GI().values())
        if method == 'get' and self.dbsummary is not None:
            return self.dbsummary
            
        tqdm.write("Getting summary from each one")
        _results = []
        if self.dbname in batch_return_dbs:
            pass
        elif self.dbname == 'protein':
            results, failed = self.edl.efetch(
                db=self.dbname,
                ids=all_GI,
                retmode="text",
                retype="gp",
                result_func=lambda x: read_efetch(x),)
            for result in results:
                self.dbsummary.update(result)
                
        elif self.dbname == 'nuccore':
            results, failed = self.edl.efetch(
                db=self.dbname,
                ids=all_GI,
                retmode="text",
                retype="gb",
                result_func=lambda x: read_efetch(x),)
            for result in results:
                self.dbsummary.update(result)
                
        elif self.dbname in single_return_dbs:
            results, failed = self.edl.esummary(
                db=self.dbname,
                ids=all_GI,
                result_func=lambda x: eread(x),
            )
            if failed:
                tqdm.write("failed retrieve summary of %s protein ID" % len(failed))
                tqdm.write("retrieve each failed GI one by one")
                _results, _failed = self.edl.esummary(
                    db=self.dbname,
                    ids=failed,
                    batch_size=1,
                    result_func=lambda x: eread(x),
                )
                tqdm.write("Following ID failed: " + "\n".join(map(str, _failed)))
            for result in results + _results:
                _dict = self.get_key_from_summary_results(result)
                self.dbsummary.update(_dict)

        elif self.dbname == "assembly":
            results, failed = self.edl.esummary(
                db="assembly", 
                ids=all_GI, 
                result_func=lambda x: parse_assembly_xml(x)
            )
            # return results is a list of defaultdict.
            for _dict in results:
                _dict = dict(_dict)
                for aid, info_dict in _dict.items():
                    info_dict["TaxId"] = info_dict["SpeciesTaxid"]
                self.dbsummary.update(_dict)
        return self.dbsummary
    
    def get_taxon(self,method='get'):
        self.check_methods(method)
        if self.dbname not in tax_convertable_dbs:
            raise IOError(f"the original ID not in '{' '.join(tax_convertable_dbs)}' ")
        dbsummary = self.get_db_summary(method=method)    
        for aid, result in dbsummary.items():
            if 'TaxId' in result:
                self.tids[aid] = int(result["TaxId"])
            elif 'taxon' in result:
                self.tids[aid] = int(result["taxon"])
        return self.tids
    
    def get_key_from_summary_results(self, retrieve_r):
        if self.dbname in ["protein", "nuccore","nucleotide"]:
            return {retrieve_r["AccessionVersion"]: retrieve_r}
    def construct_taxon_info_dict(self):
        pass

    def get_taxons_from_tid(self,tids=None):
        if tids is None:
            tids = self.tids
            # tids = self.get_taxon('get')
        # tqdm.write('get tid info using NCBITaxa')
        id2taxon = {}
        for ori_id, tid in tqdm(tids.items()):
            try:
                taxon_dict = tax2tax_info(tid)
                id2taxon[ori_id] = taxon_dict
            except ValueError as e:
                print(e)
        return id2taxon

    def update_all(self):
        # todo: update whole convertor when it add extra ID
        pass
    
    def pid2assembly(self,batch_size=50):
        if self.dbname != "protein":
            raise SyntaxError("source database must be protein")
        results, failed = self.edl.efetch(
            db="protein",
            ids=self.origin_ids,
            retype="ipg",
            retmode="text",
            batch_size= batch_size,
            result_func= lambda x: parse_ipg(x),
        )
        pid2assembly_dict = defaultdict(dict)
        for pid, nuc_info, assembly_ID in tqdm(results):
            pid2assembly_dict[pid]["assembly"] = assembly_ID
            pid2assembly_dict[pid]["nuc_info"] = nuc_info

        pid2assembly_dict.update(
            {
                _: {"assembly": "", "nuc_info": ("", "", "", "")}
                for _ in self.origin_ids
                if _ not in pid2assembly_dict
            }
        )
        return pid2assembly_dict

    def nid2assembly(self, all_GI=None):
        if self.dbname not in  ["nuccore",]:
            raise SyntaxError("source database must be nuccore")

        if all_GI is None:
            if self.GI is None:
                self.get_GI()
            all_GI = list(self.GI.values())

        results, failed = self.edl.elink(
            dbfrom="nuccore",
            db="assembly",
            ids=all_GI,
            result_func=lambda x: eread(x),
        )
        nid2assembly_dict = {}
        for r in results:
            if r["LinkSetDb"]:
                nid2assembly_dict[r["IdList"][0]] = r["LinkSetDb"][0]["Link"][0]["Id"]
        return nid2assembly_dict
        # todo: finish it

    def bp2srr(self, all_GI=None):
        # one versus multiple
        if self.dbname != "bioproject":
            raise SyntaxError("source database must be bioproject")

        if all_GI is None:
            if self.GI is None:
                self.get_GI()
            all_GI = list(self.GI.values())

        results, failed = self.edl.elink(
            dbfrom="bioproject",
            db="sra",
            ids=all_GI,
            result_func=lambda x: eread(x),
        )
        bp2srr_ids = {}
        for r in results:
            if r["LinkSetDb"]:
                bp2srr_ids[r["IdList"][0]] = [
                    _["Id"] for _ in r["LinkSetDb"][0]["Link"]
                ]
        srr_list = [srr for bp, srr_ids in bp2srr_ids.items() for srr in srr_ids]

        results
        results, failed = edl.esummary(db="sra", ids=srr_list)

        return bp2srr_ids




    
    
if __name__ == "__main__":
    pass