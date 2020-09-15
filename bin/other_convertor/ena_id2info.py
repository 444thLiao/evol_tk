# use instead !!!
# pip install enasearch
# it useless!!!!!

"""
For wgs_set
"""
import requests
import json
import enasearch

base_url = "https://www.ebi.ac.uk/ena/portal/api/links"
# url_template = f"{base_url}/sample?accession=ERS3641673&format=JSON&result=wgs_set"

def ena_sample(sample_id,result='wgs_set'):

    url = f"{base_url}/sample?accession={sample_id}&format=JSON&result={result}"
    response = requests.get(url, )
    return response

def sample2wgs(sample_id):
    response = ena_sample(sample_id,result='wgs_set')
    if response.status_code != 200:
        return ''
    else:
        out = response.text
        data = json.loads(out)
        return data[0]

def sample2assembly(sample_id):
    response = ena_sample(sample_id,result='assembly')
    if response.status_code != 200:
        return ''
    else:
        out = response.text
        data = json.loads(out)
        return data[0]



data = enasearch.search_data(
    free_text_search=True,
    query="CABMLH010000000",
    result='wgs_set',
    display='xml'
)
data = enasearch.retrieve_data(
    ids="CABMLH010000000",
    display="html")


def get_returnable_fields(result):
    url = f"https://www.ebi.ac.uk/ena/portal/api/returnFields?dataPortal=ena&format=json&result={result}"
    response = requests.get(url, )
    out = response.text
    data = json.loads(out)
    return data
returnable_fields = get_returnable_fields('wgs_set')
returnable_fields = [_['columnId'] for _ in returnable_fields]

def file_report(sample_id):
    fields = ','.join(returnable_fields)
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={sample_id}&fields={fields}&format=json&result=wgs_set"
    response = requests.get(url, )
    out = response.text
    data = json.loads(out)
    return data

returnable_fields_assembly = get_returnable_fields('assembly')
returnable_fields_assembly = [_['columnId'] for _ in returnable_fields_assembly]

def file_report_assembly(sample_id):
    fields = ','.join(returnable_fields_assembly)
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={sample_id}&fields={fields}&format=json&result=assembly"
    response = requests.get(url, )
    out = response.text
    data = json.loads(out)
    return data

from tqdm import tqdm
s2info = {}
record_failed = []
for _ in tqdm(t):
    try:
        s2info[_] = file_report(_)
    except:
        record_failed.append(_)
cmds = []
for acc,_dict in s2info.items():
    url = _dict[0]['fasta_file']
    name = url.split('/')[-1]
    cmd = f"wget {url} -O ./rawdata/extra_data/fna/{name}"
    cmds.append(cmd)
from subprocess import check_call
for cmd in cmds:
    check_call(cmd,shell=1)

s2assembly_info = {}
for _f in record_failed:
    s2assembly_info[_f] = file_report_assembly(_f)