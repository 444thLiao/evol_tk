"""
This script is mainly for retrieve infomation enough for following analysis

"""
from Bio import Entrez
from tqdm import tqdm
from Bio import Entrez, SeqIO
import io
from collections import defaultdict
import multiprocessing as mp

import requests
import io
import threading
import time

from concurrent.futures import ThreadPoolExecutor
from concurrent import futures
import pandas as pd
from bs4 import BeautifulSoup

class RequestLimiter:

    def __init__(self, min_wait=0.4):
        """The RequestLimiter class provides functionality to limit the rate at which new requests are made."""
        self.lock = threading.Lock()
        self.last_request = None
        self.min_wait = min_wait

    def wait(self):
        """The wait() function blocks until a minimum wait time from the previous invocation has passed. Thread safe."""
        with self.lock:
            # This is the first request
            if not self.last_request:
                self.last_request = time.time()
                return

            # This is not the first request
            diff = time.time() - self.last_request
            if diff < self.min_wait:
                tsleep = self.min_wait - diff
                time.sleep(tsleep)
            self.last_request = time.time()


class ResultCollector:

    def __init__(self, pbar=None):
        """The ResultCollector class provides functionality for threads to deliver their results."""
        self.pbar = pbar
        self.results = []
        self.failed = []
        self.lock = threading.Lock()

    def add_results(self, results):
        """Adds results to the collector. If a progress bar was provided, it updates the progress bar."""
        with self.lock:
            self.results += results
            if self.pbar:
                self.pbar.update(len(results))

    def add_failed(self, ids):
        """Adds failed IDs to the collector. If a progress bar was provided, it updates the progress bar."""
        with self.lock:
            self.failed += ids
            if self.pbar:
                self.pbar.update(len(ids))


class EntrezDownloader:

    def __init__(self, num_threads=30, batch_size=10, email=None, api_key=None, pbar=False):
        """The EntrezDownloader class enables parallel downloads via the NCBI Entrez interface"""
        self.baseurl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
        self.num_threads = num_threads
        self.batch_size = batch_size
        self.email = email
        self.api_key = api_key
        self.request_limiter = RequestLimiter(
            min_wait=0.37 if not api_key else 0.15)
        self.print_lock = threading.Lock()
        self.pbar = pbar

    def _efetch_batch(self, db, ids, result_collector, result_func, retmode,retype,**kwargs):
        if not retmode:
             retmode='text'
        if not retype:
            retype='gb'
        post_data = {
            'email': self.email,
            'api_key': self.api_key,
            'id': ','.join(list(map(str,ids))),
            'db': db,
            'retmode':retmode ,
            'rettype':retype 
        }

        post_data.update(kwargs)

        if self.email:
            post_data.update({'email': self.email})

        if self.api_key:
            post_data.update({'api_key': self.api_key})

        error = None
        for i in range(3):  # Retry three times
            try:
                self.request_limiter.wait()
                response = requests.post(
                    f'{self.baseurl}/efetch.cgi', post_data)
                if response.status_code == 200:
                    results = result_func(response.text)
                    result_collector.add_results(results)
                    error = None
                    break
                else:
                    error = f'[STATUS {response.status_code}] An error occurred, you may see a response text here: {response.text}'
            except Exception as e:
                error = f'[UNKNOWN ERROR] {e}'

        if error:
            result_collector.add_failed(ids)
            print(error)

    def _esummary_batch(self, db, ids, result_collector, result_func, **kwargs):
        post_data = {
            'email': self.email,
            'api_key': self.api_key,
            'id': ','.join(ids),
            'db': db,
        }

        post_data.update(kwargs)

        if self.email:
            post_data.update({'email': self.email})

        if self.api_key:
            post_data.update({'api_key': self.api_key})

        error = None
        for i in range(3):  # Retry three times
            try:
                self.request_limiter.wait()
                response = requests.post(
                    f'{self.baseurl}/esummary.cgi', post_data)
                if response.status_code == 200:
                    results = result_func(response.text)
                    result_collector.add_results(results)
                    error = None
                    break
                else:
                    error = f'[STATUS {response.status_code}] An error occurred, you may see a response text here: {response.text}'
            except Exception as e:
                error = f'[UNKNOWN ERROR] {e}'

        if error:
            result_collector.add_failed(ids)
            print(error)
    def _elink_batch(self, dbfrom,db, ids, result_collector, result_func, **kwargs):
        post_data = {
            'email': self.email,
            'api_key': self.api_key,
            'id': ids,
            'db': db,
            'dbfrom':dbfrom
        }

        post_data.update(kwargs)

        if self.email:
            post_data.update({'email': self.email})

        if self.api_key:
            post_data.update({'api_key': self.api_key})

        error = None
        for i in range(3):  # Retry three times
            try:
                self.request_limiter.wait()
                response = requests.post(
                    f'{self.baseurl}/elink.cgi', post_data)
                if response.status_code == 200:
                    results = result_func(response.text)
                    result_collector.add_results(results)
                    error = None
                    break
                else:
                    error = f'[STATUS {response.status_code}] An error occurred, you may see a response text here: {response.text}'
            except Exception as e:
                error = f'[UNKNOWN ERROR] {e}'

        if error:
            result_collector.add_failed(ids)
            print(error)
            
    def elink(self, dbfrom,db, ids, result_func=lambda x: [x], **kwargs):
        """Interface to the efetch database.
        result_func: A function to be applied to the response. Must return an iterable.
        """

        if self.pbar:
            from tqdm import tqdm
            results = ResultCollector(
                pbar=tqdm(total=len(ids), unit='records'))
        else:
            results = ResultCollector()

        executor = ThreadPoolExecutor(max_workers=self.num_threads)

        fs = []
        for start in range(0, len(ids), self.batch_size):
            num = len(ids)-start
            num = self.batch_size if num > self.batch_size else num
            f = executor.submit(self._elink_batch,
                                db=db,
                                dbfrom=dbfrom,
                                ids=ids[start:start+num],
                                result_collector=results,
                                result_func=result_func,
                                **kwargs)
            fs.append(f)

        futures.wait(fs)

        return results.results, results.failed        
    def efetch(self, db, ids, retmode,retype,result_func=lambda x: [x], batch_size=20,**kwargs):
        """Interface to the efetch database.
        result_func: A function to be applied to the response. Must return an iterable.
        """

        if self.pbar:
            from tqdm import tqdm
            results = ResultCollector(
                pbar=tqdm(total=len(ids), unit='records'))
        else:
            results = ResultCollector()

        executor = ThreadPoolExecutor(max_workers=self.num_threads)

        fs = []
        if not batch_size:
            batch_size = self.batch_size
        for start in range(0, len(ids), batch_size):
            num = len(ids)-start
            num = batch_size if num > batch_size else num
            f = executor.submit(self._efetch_batch,
                                db=db,
                                ids=ids[start:start+num],
                                result_collector=results,
                                result_func=result_func,
                                retmode=retmode,
                                retype=retype,
                                **kwargs)
            fs.append(f)

        futures.wait(fs)

        return results.results, results.failed

    def esummary(self, db, ids, result_func=lambda x: [x], **kwargs):
        """Interface to the efetch database.
        result_func: A function to be applied to the response. Must return an iterable.
        """

        if self.pbar:
            from tqdm import tqdm
            results = ResultCollector(
                pbar=tqdm(total=len(ids), unit='records'))
        else:
            results = ResultCollector()

        executor = ThreadPoolExecutor(max_workers=self.num_threads)

        fs = []
        for start in range(0, len(ids), self.batch_size):
            num = len(ids)-start
            num = self.batch_size if num > self.batch_size else num
            f = executor.submit(self._esummary_batch,
                                db=db,
                                ids=ids[start:start+num],
                                result_collector=results,
                                result_func=result_func,
                                **kwargs)
            fs.append(f)

        futures.wait(fs)

        return results.results, results.failed


edl = EntrezDownloader(
    # An email address. You might get blocked by the NCBI without specifying one.
    email='l0404th@gmail.com',
    # An API key. You can obtain one by creating an NCBI account. Speeds things up.
    api_key='ccf9847611deebe1446b9814a356f14cde08',
    num_threads=30,                       # The number of parallel requests to make
    batch_size=500,                        # The number of IDs to fetch per request
    pbar=True                             # Enables a progress bar, requires tqdm package
)

def parse_id(infile):
    id_list = []
    for row in tqdm(open(infile, 'r')):
        if row:
            id_list.append(row.split('\t')[1])
    return id_list

def batch_pid2nid(pids):
    pid2nid = {}
    try:
        nuc_ids = Entrez.read(Entrez.elink(
            dbfrom='protein', id=pids, db='nuccore'))
        pid2nid = {_['IdList'][0]: _['LinkSetDb']
                   [0]['Link'][0]['Id'] for _ in nuc_ids}
    except:
        pass
    return pid2nid

def get_elink_result(elink_results):
    _dict = {}
    for _ in elink_results:
        try:
            _dict[int(_['IdList'][0])] =  int(_['LinkSetDb'][0]['Link'][0]['Id'])
        except:
            _dict[int(_['IdList'][0])] = None
    return _dict

def batch_pid2bioproject_id(pids):
    pid2bioproject_id = {}
    bioproject_ids = Entrez.read(Entrez.elink(
            dbfrom='protein', id=pids, db='bioproject'))
    for _ in bioproject_ids:
        try:
            pid2bioproject_id[_['IdList'][0]] =  _['LinkSetDb'][0]['Link'][0]['Id']
        except:
            pid2bioproject_id[_['IdList'][0]] = None
    return pid2bioproject_id

def pid2nuc_info(pid):
    nuc_ids = Entrez.read(Entrez.elink(dbfrom='protein', id=pid, db='nuccore'))
    try:
        pid2nid = {_['IdList'][0]: _['LinkSetDb']
                   [0]['Link'][0]['Id'] for _ in nuc_ids}

        tinfo = SeqIO.read(Entrez.efetch(
            db='nuccore', id=tid, rettype='gb', retmode='text'), format='genbank')
    except:
        tid = tinfo = ''
    return tinfo


def get_info(pid):
    _cache = defaultdict(dict)

    tinfo = pid2nuc_info(pid)
    # if not tinfo:
    #     pass
    #     #failed_id.append(pid)
    # else:
    annotations = tinfo.annotations
    ref_texts = [_.title
                 for _ in annotations.get('references', [])
                 if 'Direct' not in _.title]
    _cache[pid]['reference'] = ';'.join(ref_texts)
    _cache[pid]['source'] = annotations['source']
    _cache[pid]['org'] = annotations['organism']
    _cache[pid]['keywords'] = ';'.join(annotations.get('keywords', []))
    return _cache

def batch_iter(iter, batch_size):
    # generating batch according batch_size
    n_iter = []
    batch_d = 0
    for batch_u in range(0, len(iter), batch_size):
        if batch_u != 0:
            n_iter.append(iter[batch_d:batch_u])
        batch_d = batch_u
    n_iter.append(iter[batch_d: len(iter) + 1])
    return n_iter


def parse_bioproject_xml(xml_text):
    final_project2info_dict = defaultdict(dict)
    soup = BeautifulSoup(xml_text,'xml')
    split_out = soup.find_all("DocumentSummary")
    for each_record in split_out:
        body = each_record.Project
        major_parts = {_.name:_ for _ in body.children if _.name}
        PID_part = major_parts['ProjectID']
        bioproject_id = PID_part.ArchiveID['accession']
        bioproject_gid = int(PID_part.ArchiveID['id'])
        PID_des = major_parts['ProjectDescr']
        des_text = PID_des.Description
        des_title = PID_des.Title
        PID_pubmed_id = PID_des.Publication['id']
        list_biosample = PID_des.find_all('LocusTagPrefix')
        biosample_text = ';'.join([_['biosample_id'] for _ in list_biosample])
        
        PID_type = major_parts['ProjectType']
        biological_properties = PID_type.find('BiologicalProperties')
        if not biological_properties:
            env_data = biological_properties.find('Environment')
        
        final_project2info_dict[bioproject_id]['GI'] = bioproject_gid
        final_project2info_dict[bioproject_id]['description'] = des_text
        final_project2info_dict[bioproject_id]['title'] = des_title
        final_project2info_dict[bioproject_id]['pubmed GI'] = PID_pubmed_id
        final_project2info_dict[bioproject_id]['relative biosample'] = biosample_text
        final_project2info_dict[bioproject_id]['number of biosamples'] = len(list_biosample)
        for _data in env_data:
            if _data !='\n':
                final_project2info_dict[bioproject_id][_data.name] = _data.text
    return final_project2info_dict

def main(infile):
    id_list = parse_id(infile)
    id_list = list(set(id_list))
    pid2info_dict = defaultdict(dict)
    
    
    tqdm.write('get pid summary from each one')
    results, failed = edl.esummary(db='protein',
                                   ids=id_list,
                                   result_func=lambda x: Entrez.read(
                                       io.StringIO(x)))
    gi2pid = {}
    for result in results:
        aid = result['AccessionVersion']
        gi = result['Gi'].real
        taxid = result['TaxId'].real
        pid2info_dict[aid]['GI'] = gi
        pid2info_dict[aid]['taxid'] = taxid
        gi2pid[gi] = aid
        
    tqdm.write('retrieving nuccore info')
    prot_results, prot_failed = edl.efetch(db='protein',
                                   ids=list(map(str,id_list)),
                                   result_func=lambda x: list(SeqIO.parse(io.StringIO(x),format='genbank'))
                                   )
    
    for prot_t in prot_results:
        aid = prot_t.id
        if aid not in id_list:
            print('error ', aid)
        annotations = prot_t.annotations
        
        ref_texts = [_.title
                 for _ in annotations.get('references', [])
                 if 'Direct' not in _.title]
        for idx,ref_t in enumerate(ref_texts):
            pid2info_dict[aid]['reference_'+str(int(idx)+1)] = ref_t
        pid2info_dict[aid]['nuccore id'] = annotations['db_source'].split(' ')[-1]
        pid2info_dict[aid]['source'] = annotations['source']
        pid2info_dict[aid]['org'] = annotations['organism']
        pid2info_dict[aid]['keywords'] = ';'.join(annotations.get('keywords', []))
        pid2info_dict[aid]['comments'] = annotations.get('comment','')
        
        pid2info_dict[aid].update(dict([_.split(':') for _ in prot_t.dbxrefs]))
    
    
    tqdm.write('processing pid to bioproject id and biosample')
    results, failed = edl.efetch(db='bioproject',
                                ids=[d.get('BioProject','') for d in pid2info_dict.values() if 'BioProject' in d][:5],
                                retmode='xml',
                                retype='xml',
                                #result_func=lambda x: Entrez.read(io.StringIO(x),validate=False)
                                )
    
    
    
    pid2bioproject_dict = get_elink_result(results)
    for pro_gi,bioporject_gi in pid2bioproject_dict.items():
        pid2info_dict[gi2pid[pro_gi]]['bioproject gi'] = bioporject_gi
        
        
    final_df = pd.DataFrame.from_dict(pid2info_dict,orient='index')
    return final_df
    
    
    
    
    
    tqdm.write('get pid to nid')
    results, failed = edl.elink(dbfrom='protein',
                                ids=id_list,
                                db='nuccore',
                                result_func=lambda x: Entrez.read(io.StringIO(x))
                                )
    pid2nid_dict = get_elink_result(results)
    for pro_gi,nuc_gi in pid2nid_dict.items():
        pid2info_dict[gi2pid[int(pro_gi)]]['nuccore gi'] = int(nuc_gi)
        
        
    tqdm.write('get nid summary from each one')
    
    



    nid_list = [pid2info_dict[d]['nuccore gi'] for d in pid2info_dict.keys()]

    
    for nuccore_t in nuccore_results:
        annotations = nuccore_t.annotations
        ref_texts = [_.title
                 for _ in annotations.get('references', [])
                 if 'Direct' not in _.title]
        _cache[pid]['reference'] = ';'.join(ref_texts)
        _cache[pid]['nuccore id'] = nuccore_t.id
        _cache[pid]['source'] = annotations['source']
        _cache[pid]['org'] = annotations['organism']
        _cache[pid]['keywords'] = ';'.join(annotations.get('keywords', []))
        pid2info_dict[gi2pid[pro_gi]]['bioproject gi']
