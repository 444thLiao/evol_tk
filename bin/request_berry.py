import requests
from tqdm import tqdm
import json,io
import urllib
from os.path import *
import warnings

def download_file(url,local_filename):
    # NOTE the stream=True parameter below
    with requests.get(url, stream=True,
                      verify=False) as r:
        r.raise_for_status()
        total_size = int(r.headers.get('content-length', 0))
        with open(local_filename, 'wb') as f:
            for chunk in tqdm(r.iter_content(chunk_size=1024*1024),
                              total=int(total_size/(1024*1024))): 
                if chunk: # filter out keep-alive new chunks
                    f.write(chunk)
                    f.flush()
    # return local_filename
warnings.filterwarnings('ignore')
post_data = {'func':'get_list',
               'ssid':'02D9JBk',
               'fid':'02D9JBk',
               'path':'/191101_X515_0555_AH2MVCCCX2',
               'sort':'natural',
               'dir':'ASC'}
r = requests.post('https://ss.berrygenomics.com/share.cgi',
              post_data,
              verify=False)
data = json.load(io.StringIO(r.text))
total = data['total']
start = 0
done_count = 0
while done_count <= total:
    for f_dict in tqdm(data['datas']):
        filename = f_dict['filename']
        _p = dict(ssid=post_data['ssid'],
                fid=post_data['fid'],
                path=post_data['path'],
                filename=filename,
                openfolder='forcedownload',
                ep='',
                start=start)
        url_r = "https://ss.berrygenomics.com/share.cgi?" + urllib.parse.urlencode(_p)
        if exists(filename):
            if int(f_dict['filesize']) == getsize(filename):
                continue
        tqdm.write('downloading %s' % filename)
        download_file(url_r,filename)
        done_count +=1
    start+=50
    post_data.update({'start':start})
    r = requests.post('https://ss.berrygenomics.com/share.cgi',
              post_data,
              verify=False)
    data = json.load(io.StringIO(r.text))
    
    
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        raise Exception("UNKNOWN parameters... first is url, second is local_filename")
    download_file(sys.argv[1],
                  sys.argv[2])
    
# url = "https://ss.berrygenomics.com/share.cgi?ssid=02D9JBk&fid=02D9JBk&path=%2F&filename=191101_X515_0555_AH2MVCCCX2&openfolder=forcedownload&ep="
# download_file(url,'download_2019-11-13_15-00-51.zip')