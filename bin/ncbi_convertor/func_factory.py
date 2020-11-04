from bin.ncbi_convertor.pid2GI import main
import io
from collections import defaultdict

from Bio import Entrez, SeqIO
from tqdm import tqdm

from bin.ncbi_convertor.toolkit import edl, access_intermedia, tax2tax_info, parse_ipg
from global_search.thirty_party.metadata_parser import parse_assembly_xml

tax_convertable_dbs = ['protein', 'assembly', 'nuccore']
batch_return_dbs = []
single_return_dbs = ['protein', 'nuccore']

edl.num_threads = 5
class NCBI_convertor():
    "prototype of function which convert ID to GI"

    def __init__(self, IDs, db):
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

    def load_from_file(self,infile):
        # todo:
        pass
    def print_available_dbs(self):
        print("")

    def get_biosample(self):
        pass

    def get_bioproject(self):
        pass

    def check_cache(self, suffix, redo):
        _cache = access_intermedia(self.origin_ids,
                                   suffix=suffix, redo=redo)
        if _cache is not None:
            id2other_IDs = _cache
            if "2GI" in suffix:
                self.GI = id2other_IDs
            elif '' in suffix:
                pass
        else:
            print("No cache found.")

    def get_GI(self, num_retry=5):
        if self.GI is not None:
            return
        tqdm.write("Get GI......")
        results, failed = edl.esearch(db=self.dbname,
                                      ids=self.origin_ids,
                                      result_func=lambda x: Entrez.read(
                                          io.StringIO(x))['IdList'],
                                      batch_size=1
                                      )
        # manual retry the retrieval in case of network problems
        _results_dict = {}
        _count = 0
        while failed:
            failed_id_list = failed
            _results, failed = edl.esearch(db=self.dbname,
                                           ids=failed_id_list,
                                           result_func=lambda x: Entrez.read(
                                               io.StringIO(x))['IdList'],
                                           batch_size=1
                                           )
            _results_dict.update(dict(_results))
            _count += 1
            if _count >= num_retry:
                break

        tqdm.write('still %s failed IDs, be careful.....' % len(failed))
        # for edl.esearch, it will auto **zip** searched term and its result.
        id2gi = dict(results)
        id2gi.update(_results_dict)
        id2gi = {pid: id2gi.get(pid, '') for pid in self.origin_ids}
        # stodge the result into intermedia file for second access.
        process_name = f"{self.dbname}2GI"
        access_intermedia(id2gi, suffix=process_name)

        self.GI = id2gi
        self.failed_ids = failed

    def get_db_summary(self):
        if self.GI is None:
            self.get_GI()

        all_GI = list(set(self.GI.values()))
        tqdm.write('Getting summary from each one')
        _results = []
        if self.dbname in batch_return_dbs:
            pass

        elif self.dbname in single_return_dbs:
            results, failed = edl.esummary(db=self.dbname,
                                           ids=all_GI,
                                           result_func=lambda x: Entrez.read(
                                               io.StringIO(x))
                                           )

            if failed:
                tqdm.write("failed retrieve summary of %s protein ID" % len(failed))
                tqdm.write("retrieve each failed GI one by one")
                _results, _failed = edl.esummary(db=self.dbname,
                                                 ids=failed,
                                                 batch_size=1,
                                                 result_func=lambda x: Entrez.read(
                                                     io.StringIO(x))
                                                 )
                tqdm.write("Following ID failed: " + '\n'.join(map(str, _failed)))
            for result in results + _results:
                _dict = self.get_key_from_summary_results(result)
                self.dbsummary.update(_dict)

        elif self.dbname == 'assembly':
            results, failed = edl.esummary(db='assembly',
                                           ids=all_GI,
                                           result_func=lambda x: parse_assembly_xml(x))
            # return results is a list of defaultdict.
            for _dict in results:
                _dict = dict(_dict)
                for aid, info_dict in _dict.items():
                    info_dict['TaxId'] = info_dict['SpeciesTaxid']
                self.dbsummary.update(_dict)

    def get_taxon(self, ):
        if self.dbname not in tax_convertable_dbs:
            raise IOError(f"the original ID not in ")
        if not self.dbsummary:
            self.get_db_summary()
        for aid, result in self.dbsummary.items():
            self.tids[aid] = int(result['TaxId'])

    def get_key_from_summary_results(self, retrieve_r):
        if self.dbname in ['protein', 'nuccore']:
            return {retrieve_r['AccessionVersion']: retrieve_r}

    def construct_taxon_info_dict(self):
        pass

    def get_taxons_from_tid(self):
        if not self.tids:
            self.get_taxon()
        id2taxon = {}
        for ori_id, tid in self.tids.items():
            taxon_dict = tax2tax_info(tid)
            id2taxon[ori_id] = taxon_dict
        return id2taxon

    def update_all(self):
        # todo: update whole convertor when it add extra ID
        pass

    def pid2assembly(self):
        if self.dbname != 'protein':
            raise SyntaxError("source database must be protein")
        self.get_GI()
        all_GI = list(self.GI.values())
        results, failed = edl.efetch(db='protein',
                                     ids=all_GI,
                                     retmode='ipg',
                                     retype='xml',
                                     result_func=lambda x: parse_ipg(x))
        pid2assembly_dict = defaultdict(dict)
        for pid, nuc_info, assembly_ID in tqdm(results):
            pid2assembly_dict[pid]['assembly'] = assembly_ID
            pid2assembly_dict[pid]['nuc_info'] = nuc_info

        pid2assembly_dict.update({_: {'assembly': '',
                                      'nuc_info': ('', '', '', '')
                                      }
                                  for _ in self.origin_ids
                                  if _ not in pid2assembly_dict})
        return pid2assembly_dict

    def nid2assembly(self):
        if self.dbname != 'nuccore':
            raise SyntaxError("source database must be nuccore")
        self.get_GI()
        all_GI = list(self.GI.values())
        # results, failed = edl.efetch(db='nuccore',
        #                              ids=all_GI,
        #                              retmode='text',
        #                              retype='gb',
        #                              result_func=lambda x: list(SeqIO.parse(io.StringIO(x), format='genbank')))
        
        
        results, failed = edl.elink(dbfrom='nuccore',
                                    db='assembly',
                                ids=all_GI,
                                batch_size=1,
                                result_func=lambda x: Entrez.read(
                                                     io.StringIO(x)))
        nid2assembly_dict = {}
        for r in results:
            if r['LinkSetDb']:
                nid2assembly_dict[r['IdList'][0]] = r['LinkSetDb'][0]['Link'][0]['Id']
        
        return nid2assembly_dict
        # todo: finish it

    def transform_into_other_db(self, another_db):
        pass

if __name__ == "__main__":
    tab = pd.read_csv("/home-user/jjtao/Rhizobiales/FLnif-query/gene/query_result/nifH_custom.blast",sep='\t',header=None)
    all_nuc_id = tab[0]
    all_nuc_id = list(set(all_nuc_id))
    nc = NCBI_convertor(all_nuc_id,'nuccore')
    nc.get_GI()
    all_GI = list(nc.GI.values())
    
    results, failed = edl.elink(dbfrom='nuccore',
                                db='assembly',
                            ids=all_GI,
                            batch_size=1,
                            result_func=lambda x: Entrez.read(
                                                    io.StringIO(x)))
    have_assembly_nuc = {}
    for r in results:
        if r['LinkSetDb']:
            have_assembly_nuc[r['IdList'][0]] = r['LinkSetDb'][0]['Link'][0]['Id']
    
    all_genome_id = [v for k,v in have_assembly_nuc.items() if v]
    gid2assembly_info, bp2info, bs2info = genomeID2Bio(all_genome_id)
    
    ginfo_df = pd.DataFrame.from_dict(gid2assembly_info, orient='index')
    # ginfo_df.index = ginfo_df.iloc[:,0]
    bp_df = pd.DataFrame.from_dict(bp2info, orient='index')
    bs_df = pd.DataFrame.from_dict(bs2info, orient='index')
    _df1 = bp_df.reindex(ginfo_df.loc[:, 'BioprojectAccn'])
    _df1.index = ginfo_df.index
    _df2 = bs_df.reindex(ginfo_df.loc[:, 'BioSampleAccn'])
    _df2.index = ginfo_df.index
    full_df = pd.concat([ginfo_df,
                            _df1,
                            _df2], axis=1)
    full_df = full_df.applymap(lambda x: x.replace('\n', ' ')
    if isinstance(x, str) else x)
    full_df = full_df.drop(['GI', 'relative biosample'], axis=1)