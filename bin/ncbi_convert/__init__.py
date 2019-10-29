"""
Entry for using ncbi convertor.

"""
from global_search.thirty_party.EntrezDownloader import EntrezDownloader
from global_search.classification_script import _classificated
from global_search.thirty_party.metadata_parser import *
from tqdm import tqdm
import hashlib
from os.path import exists,dirname,join
import json
import os


tmp_dir = './.tmp_getINFO'
taxons = ['superkingdom', 'phylum', 'class',
          'order', 'family', 'genus', 'species']
import hashlib
def shash(astr):
    return hashlib.md5(astr.encode()).hexdigest()

def access_intermedia(obj,suffix=''):
    """
    It is mainly stodge intermediate result from previous results
    or provide a way to access it.
    Normally it should archieve into current directory with md5 hashed file name.
    
    obj is necessary, normally it is a list of IDs/dictionary which also needed to genereate md5 hashed file name.
    ofile is optional depend on what you want
    """
    nameofid = list(set(obj))
    _md5 = str(shash(';'.join(list(sorted(nameofid)))  ))
    ofile = join(tmp_dir,_md5) +'_' +suffix
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile),exist_ok=1)
    
    if exists(ofile):
        load_obj = json.load(open(ofile,'r'))
        if isinstance(load_obj,dict):
            return load_obj
    else:
        if isinstance(obj,list):
            # just used to validate run it yet?
            pass
        elif isinstance(obj,dict):
            json.dump(obj,open(ofile,'w'))
        else:
            raise Exception('unexpected data type')

def parse_id(infile, columns=0):
    """
    1. format protein accession ID instead of using wrong one
    2. retrieving information provided inside the infile
    """
    id_list = []
    id2info = {}
    header_info = ''
    for row in tqdm(open(infile, 'r')):
        if row.startswith('#'):
            # header row
            header_info = row.strip('\n').strip('#').split('\t')
        if row:
            id = row.split('\t')[columns].strip().strip('\n')
            if '||' in id:
                id = id.split('||')[-1]
            id_list.append(id)
            if not header_info:
                id2info[id] = ';'.join(
                row.split('\t')[columns+1:]).strip('\n')
            else:
                id2info[id] = dict(zip(header_info[1:],
                row.strip('\n').split('\t')[columns+1:]))
    return id_list, id2info



edl = EntrezDownloader(
    # An email address. You might get blocked by the NCBI without specifying one.
    email='l0404th@gmail.com',
    # An API key. You can obtain one by creating an NCBI account. Speeds things up.
    api_key='ccf9847611deebe1446b9814a356f14cde08',
    num_threads=30,   # The number of parallel requests to make
    # The number of IDs to fetch per request
    batch_size=50,
    pbar=True  # Enables a progress bar, requires tqdm package
)