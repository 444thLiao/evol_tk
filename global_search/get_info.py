"""
This script is mainly for retrieve infomation enough for following analysis

"""
from tqdm import tqdm


def parse_id(infile):
    id_list = []
    for row in tqdm(open(infile,'r')):
        if row:
            id_list.append(row.split('\t'))
    return id_list

def main(infile):
    id_list = parse_id(infile)
    
    
