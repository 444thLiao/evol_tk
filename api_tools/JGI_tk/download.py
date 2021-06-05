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
  'Cookie': """externallySequenced="true"; showRestricted="true"; activateHighlights="false"; showOnlyPublished="false"; showSuperseded="false"; showAll="false"; searchWebExecuted="false"; JSESSIONID=BFD9B5A8922A32C86DDC96B0D34E7E2F; currentUserId="l0404th@gmail.com"; query="2721755844|2718217642|2708742552|2654588083|2654587960|2619618950|2558309099|2545555825|2527291509|2527291500|2524023104|2513237068|2513237066|2264867229|2264867070|2264867067|2263082000"; genomeSearchParams="?core=genome&query=2721755844%7C2718217642%7C2708742552%7C2654588083%7C2654587960%7C2619618950%7C2558309099%7C2545555825%7C2527291509%7C2527291500%7C2524023104%7C2513237068%7C2513237066%7C2264867229%7C2264867070%7C2264867067%7C2263082000'&searchIn=Anything&searchType=Keyword&showAll=false&externallySequenced=true&sortBy=displayNameStr&showRestricted=true&showOnlyPublished=false&showSuperseded=true&sortOrder=asc&rawQuery=false&showFungalOnly=false&activateHighlights=false&programName=all&programYear=all&superkingdom=--any--&status=--any--&scientificProgram=--any--&productCategory=--any--"; searchGenomeExecuted="false"; savedSearchNames=[{"lastUpdated":"2021-06-04 18:33:02","name":"tmp","fullName":"tmp"}]; _ga=GA1.2.1287611267.1622814897; _gid=GA1.2.173042913.1622814897; __utmc=249693981; __utmz=249693981.1622816314.1.1.utmcsr=(direct)|utmccn=(direct)|utmcmd=(none); jgi_session=/api/sessions/6b71279de53c0b18934e8cbe4f1fdcf4; __utma=249693981.1287611267.1622814897.1622816314.1622856202.2; jgi_return=https://img.jgi.doe.gov//cgi-bin/mer/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=2545555825; __utmb=249693981.14.9.1622858562261""",
}

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