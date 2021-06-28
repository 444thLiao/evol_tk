"""
Provide utilies for retrieving information relative to MG rast 


"""
import json
import requests
import io
import pandas as pd
def json2df(infos):
    new_infos = []
    for info in infos:
        if 'attributes' in info:
            attr = info.pop('attributes')
            for k,v in attr.items():
                info[f"attr_{k}"] = v
            new_infos.append(info)
    return pd.DataFrame.from_dict(new_infos)
def project_info(project_id):
    API_URL= f"https://api.mg-rast.org/project/{project_id}?verbosity=full"
    response = requests.get(API_URL)    
    json_text = response.text
    project_info = json.loads(json_text)
    samples = dict(project_info['samples'])
    metagenomes = project_info['metagenomes']
    info_df = json2df(metagenomes)
    metagenome_ids = list(info_df['metagenome_id'])
    return info_df


def get_MG_rast_files(mgid):
    API_URL= f"https://api.mg-rast.org/1/download/{mgid}"
    response = requests.get(API_URL)    
    json_text = response.text
    download_info = json.loads(json_text)['data']
    seq_list = [_ for _ in download_info if _.get('data_type','')=='sequence']
    fastq = [_ for _ in seq_list if _.get('file_format','') == 'fastq']
    
    links = []
    if fastq:
        links.append((fastq[0]['url'],fastq[0]['file_name']))
    return links

if __name__ =='__main__':
    pid = "mgp3359"
    info_df = project_info(pid)
    links = []
    for mgid in list(info_df['metagenome_id']):
        links.extend(get_MG_rast_files(mgid))
    
    
    
    