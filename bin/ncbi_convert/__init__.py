"""
Entry for using ncbi convertor.

"""
from global_search.thirty_party.EntrezDownloader import EntrezDownloader
from global_search.classification_script import _classificated
from global_search.thirty_party.metadata_parser import *
from tqdm import tqdm
import hashlib
from os.path import exists,dirname,join,realpath,expanduser
import json
import os


tmp_dir = '~/.tmp_getINFO'
tmp_dir = expanduser(tmp_dir)
taxons = ['superkingdom', 'phylum', 'class',
          'order', 'family', 'genus', 'species']
import hashlib
def shash(astr):
    return hashlib.md5(astr.encode()).hexdigest()

def access_intermedia(obj,suffix='',redo=False):
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
    
    if exists(ofile) and not redo:
        load_obj = json.load(open(ofile,'r'))
        if isinstance(load_obj,dict):
            tqdm.write('Dectect same cache, use it instead of run it again.')
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
            continue
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

def unpack_gb(prot_t):
    """
    For standarizing unpack process from genbank
    """
    if not prot_t:
        return {}
    _cache_dict = {}
    annotations = prot_t.annotations
    ref_texts = [_
                         for _ in annotations.get('references', [])
                         if 'Direct' not in _.title and _.title]
    seq = str(prot_t.seq)
    org = annotations.get('organism', '')
    source = annotations.get('source', '')
    db_ = dict([_.split(':') for _ in prot_t.dbxrefs if ':' in _])
    biop = db_.get('BioProject', '')
    bios = db_.get('BioSample', '')
    nuccore_ID = annotations.get('db_source', '').split(' ')[-1]
    kw = ';'.join(annotations.get('keywords', []))
    comment = annotations.get('comment', '')
    for idx in range(10):
        if idx < len(ref_texts):
            ref_t = ref_texts[idx]
            _cache_dict[f'reference title {idx+1}'] = ref_t.title
            _cache_dict[f'reference journal {idx+1}'] = ref_t.journal
            _cache_dict[f'reference authors {idx+1}'] = ref_t.authors
    _cache_dict['organism'] = org
    _cache_dict['sequence'] = seq
    _cache_dict['source'] = source
    _cache_dict['bioproject'] = biop
    _cache_dict['biosample'] = bios
    _cache_dict['nuccore ID'] = nuccore_ID
    _cache_dict['keyword'] = kw
    _cache_dict['comment'] = comment
    return _cache_dict
    
    
    
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