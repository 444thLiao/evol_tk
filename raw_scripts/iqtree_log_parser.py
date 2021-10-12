import os
import re


def log_parser(log):
    full_text = open(log).read()
    seqs = [_ for _ in full_text.split('\n') if _.startswith('Alignment has')][0]
    num_seqs,sites = seqs.split(' ')[2], seqs.split(' ')[5]
    
    RAM_estimated = [_ for _ in full_text.split('\n') if _.startswith('NOTE:')]
    RAM_estimated = ''.join(RAM_estimated[0].split(' ')[1:3])
    
    




