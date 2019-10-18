"""
This script is mainly for parse diamond output and get result ID from big database.

"""

from glob import glob 
from tqdm import tqdm
import click
import os
from os.path import join,basename,dirname,exists
from collections import defaultdict
header = 'qseqid sseqid salltitles pident length evalue bitscore'.split(' ')

def parse_multi(in_dir,suffix='faa'):
    id2name = {}
    for in_file in glob(join(in_dir,f'*.{suffix}')):
        name = basename(in_file).replace(f'.{suffix}','')
        for _ in open(in_file):
            if _.startswith('>'):
                id2name[_.strip('>').split(' ')[0].strip('\n')] = name
    return id2name           

def main(in_file,id2name = {}):
    num_line = int(os.popen(f'wc -l {in_file}').read().split(' ')[0].strip())
    sid2qids = defaultdict(set)    
    for row in tqdm(open(in_file,'r'),total=num_line):
        h2info = dict(zip(header,row.split('\t')))
        qid = h2info['qseqid']
        salltitles = h2info['salltitles']
        IDs = [_.split(' ')[0] for _ in salltitles.split('<>')]
        if not id2name:
            name = qid
        else:
            name = id2name[qid]
        for id in IDs:
            sid2qids[id].add(name)
    return sid2qids

@click.command()
@click.option('-i','infile')
@click.option('-o','ofile')
@click.option('-in_dir','--parse_indr','extra_indir',required=False)
def cli(infile,ofile,extra_indir):
    if (not exists(dirname(ofile))) and ('/' in ofile):
        os.makedirs(dirname(ofile))
    if extra_indir:
        id2name=parse_multi(extra_indir)
    else:
        id2name = {}
    sid2qids = main(infile,id2name)
    with open(ofile,'w') as f1:
        for id,vals in sid2qids.items():
            vals_str = '\t'.join(sorted(list(vals)))
            print(f"{id}\t{vals_str}",file=f1)
    

if __name__ == "__main__":
    cli()