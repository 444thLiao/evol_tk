"""
raw script for JGI download

"""

from collections import defaultdict
import requests
from bs4 import BeautifulSoup
from tqdm import tqdm
import multiprocessing as mp
from os.path import *

headers = {
  'Cookie': """ """,
}

# def get_portal_name(img_id):
#     # some img_id is different to its protal name
    
#     # the web content is dynamic generated. hard to retrieve
#     url = f"https://genome.jgi.doe.gov/portal/?core=genome&query={img_id}"
#     response = requests.get(url,headers=headers)    
#     xml_text = response.text    
#     soup = BeautifulSoup(xml_text)
#     doms = soup.find_all('tbody')
#     for b in doms[0].find_all('b'):
#         if b.text == 'Project: ':
#             break

img_ids = ['IMG_2721755844',
 'IMG_2718217642',
 'IMG_2708742552',
 'IMG_2654588083',
 'IMG_2654587960',
 'IMG_2619618950',
 'IMG_2558309099',
 'IMG_2545555825',
 'IMG_2527291509',
 'IMG_2527291500',
 'IMG_2524023104',
 'IMG_2513237068',
 'IMG_2513237066',
 'IMG_2264867229',
 'IMG_2264867070',
 'IMG_2264867067',
 'IMG_2263082000']

img_id2url = {}
img_id2raws = defaultdict(list)
for img_id in img_ids:
    url = f"https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism={img_id}&organizedByFileType=false"
    response = requests.get(url,headers=headers)    
    xml_text = response.text
    soup = BeautifulSoup(xml_text)
    for d in soup.find_all('file'):
        if d.attrs['filename'] == img_id.split('_')[-1]+'.tar.gz':
            img_id2url[img_id] = "https://genome.jgi.doe.gov" + d.attrs['url']
        else:
            img_id2raws[img_id].append(d)

def get_id(img_id,url):
    if exists(f"./{img_id}.tar.gz"):
        return 
    r = requests.get(url,headers=headers)
    print(img_id,r.status_code)
    with open(f"./{img_id}.tar.gz",'wb') as f:
        f.write(r.content)    

def run(args):
    get_id(*args)


params = [(img_id,url) for img_id,url in img_id2url.items()]
with mp.Pool(processes=20) as tp:
    r = list(tqdm(tp.imap(run, params), total=len(params)))


remaining_ids = [_ for _ in img_ids if _ not in img_id2url]
missing_ids = {"IMG_2619618950":"AigarcMDMJNZ1K18_FD",
"IMG_2527291509":"candivSAAA471B22_FD",
"IMG_2527291500":"MargroSAAA007O23_FD",
"IMG_2264867229":"candiva000106J15_FD",
"IMG_2264867070":"candivSAAA011G17_FD",
"IMG_2264867067":"candivSAAA011E11_FD",}

img_id2url = {}
for img_id,org in missing_ids.items():
    url = f"https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism={org}&organizedByFileType=false"
    response = requests.get(url,headers=headers)    
    xml_text = response.text
    soup = BeautifulSoup(xml_text)
    
    for d in soup.find_all('file'):
        if d.attrs['filename'] == img_id.split('_')[-1]+'.tar.gz':
            img_id2url[img_id] = "https://genome.jgi.doe.gov" + d.attrs['url']

params = [(img_id,url) for img_id,url in img_id2url.items()]
with mp.Pool(processes=20) as tp:
    r = list(tqdm(tp.imap(run, params), total=len(params)))