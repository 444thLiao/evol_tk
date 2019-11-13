import requests
from tqdm import tqdm

def download_file(url,local_filename):
    # NOTE the stream=True parameter below
    with requests.get(url, stream=True,verify=False) as r:
        r.raise_for_status()
        total_size = int(r.headers.get('content-length', 0))
        with open(local_filename, 'wb') as f:
            for chunk in tqdm(r.iter_content(chunk_size=1024*1024),total=total_size/(1024*1024)): 
                if chunk: # filter out keep-alive new chunks
                    f.write(chunk)
                    # f.flush()
    # return local_filename
    

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        raise Exception("UNKNOWN parameters... first is url, second is local_filename")
    download_file(sys.argv[1],
                  sys.argv[2])
    
# url = "https://ss.berrygenomics.com/share.cgi?ssid=02D9JBk&fid=02D9JBk&path=%2F&filename=191101_X515_0555_AH2MVCCCX2&openfolder=forcedownload&ep="
# download_file(url,'download_2019-11-13_15-00-51.zip')