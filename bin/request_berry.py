"""
This script is mainly for download genome from berry genomes company.
"""
import io
import json
import urllib
import warnings
from os.path import *

import click
import requests
from tqdm import tqdm


def download_file(url, local_filename):
    # NOTE the stream=True parameter below
    with requests.get(url, stream=True,
                      verify=False) as r:
        r.raise_for_status()
        total_size = int(r.headers.get('content-length', 0))
        with open(local_filename, 'wb') as f:
            for chunk in tqdm(r.iter_content(chunk_size=1024 * 1024),
                              total=int(total_size / (1024 * 1024))):
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)
                    f.flush()
    # return local_filename


def main(ssid, path, dry_run):
    warnings.filterwarnings('ignore')
    post_data = {'func': 'get_list',
                 'ssid': ssid,
                 'fid': ssid,
                 'path': f'/{path}',
                 'sort': 'natural',
                 'dir': 'ASC'}
    r = requests.post('https://ss.berrygenomics.com/share.cgi',
                      post_data,
                      verify=False)
    d = []
    data = json.load(io.StringIO(r.text))
    total = data['total']
    start = 0
    done_count = 0
    tqdm.write("retrieving all download info, please wait a moment...")
    while done_count < total:
        for f_dict in data['datas']:
            filename = f_dict['filename']
            _p = dict(ssid=post_data['ssid'],
                      fid=post_data['fid'],
                      path=post_data['path'],
                      filename=filename,
                      openfolder='forcedownload',
                      ep='',
                      start=start)
            url_r = "https://ss.berrygenomics.com/share.cgi?" + urllib.parse.urlencode(_p)

            if exists(filename) and int(f_dict['filesize']) == getsize(filename):
                done_count += 1
                continue
            d.append((url_r, filename))
            done_count += 1
        start += 50
        post_data.update({'start': start})
        r = requests.post('https://ss.berrygenomics.com/share.cgi',
                          post_data,
                          verify=False)
        data = json.load(io.StringIO(r.text))
    # start downloading
    for url_r, filename in tqdm(d):
        if dry_run:
            continue
        download_file(url_r, filename)
        tqdm.write('downloading %s' % filename)
    tqdm.write("finishing all......")


@click.command()
@click.argument("ssid")
@click.argument("path")
@click.option("-d", "--dry_run", "dry_run", default=False, is_flag=True)
def cli(ssid, path, dry_run):
    main(ssid, path, dry_run)


if __name__ == "__main__":
    cli()

# url = "https://ss.berrygenomics.com/share.cgi?ssid=02D9JBk&fid=02D9JBk&path=%2F&filename=191101_X515_0555_AH2MVCCCX2&openfolder=forcedownload&ep="
# download_file(url,'download_2019-11-13_15-00-51.zip')
# python request_berry.py 01VCdO3 200109_X499_0522_AH3K2VCCX2 --dry-run
