"""
This script is mainly for downloading genome from berry genomes sequencing company.
"""

import io
import json
import os
import urllib
import warnings
from os.path import *

import click
import requests
from tqdm import tqdm


def process_path(path):
    if not '/' in path:
        path = './' + path
    if path.startswith('~'):
        path = expanduser(path)
    if path.startswith('.'):
        path = abspath(path)
    return path


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


def main(ssid, odir, dry_run):
    odir = process_path(odir)
    if not exists(odir):
        os.makedirs(odir)
    warnings.filterwarnings('ignore')
    post_data = {'func': 'get_list',
                 'ssid': ssid,
                 'fid': ssid,
                 # 'path': f'/{path}',
                 'sort': 'natural',
                 'dir': 'ASC',
                 'ep': ''}
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
            if filename == '002':
                done_count += 1
                continue
            _p = dict(ssid=post_data['ssid'],
                      fid=post_data['fid'],
                      path='/',
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
        download_file(url_r, join(odir, filename))
        tqdm.write('downloading %s' % filename)
    tqdm.write("finishing all......")


@click.command()
@click.argument("ssid")
@click.argument("odir", default='./')
@click.option("-d", "--dry_run", "dry_run", default=False, is_flag=True)
def cli(ssid, odir, dry_run):
    main(ssid, odir, dry_run)


if __name__ == "__main__":
    cli()

# url = "https://ss.berrygenomics.com/share.cgi?ssid=02D9JBk&fid=02D9JBk&path=%2F&filename=191101_X515_0555_AH2MVCCCX2&openfolder=forcedownload&ep="
# download_file(url,'download_2019-11-13_15-00-51.zip')
# python request_berry.py 01VCdO3 200109_X499_0522_AH3K2VCCX2 --dry-run
