import pandas as pd
from os.path import join,exists,basename,dirname,isdir
from glob import glob
import os

reference_file = './ref_id'
target_dir = './genbank'

id_list = open(reference_file).read().split('\n')
for each_id in id_list:
    if not glob(join(target_dir,'*',each_id)):
        print(each_id)


target_dir = './genome_protein_files'

id_list = open(reference_file).read().split('\n')
for each_id in id_list:
    if not glob(join(target_dir,each_id+'*')):
        print(each_id)