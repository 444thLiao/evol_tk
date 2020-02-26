#################################################################################
####  Extract protein sequences from output of orthofinder
####  This script is mainly focus on the number of genomes found at single orthogroup
####
#################################################################################
import sys
import os
from os.path import dirname

sys.path.insert(0, dirname(dirname(__file__)))
import pandas as pd
from Bio import SeqIO
from os.path import join, dirname
import click
from toolkit.utils import get_dict, get_protein, get_single_copy, get_summary_statistic
from tqdm import tqdm
import multiprocessing as mp

def extract2file(args):
    sub_data,output_dir,og,species_path_temp,id2spe,id2seq = args
    with open(join(output_dir, og + '.faa'), 'w') as f1:
        seqs = []
        for speid, seq_id in sub_data.loc[og, :].items():
            if speid != seq_id.split('_')[0] and seq_id != 'nan':
                speid = seq_id.split('_')[0]
            spe_file = species_path_temp.format(speid=speid)
            if not os.path.exists(spe_file):
                _cache = spe_file.rpartition('WorkingDirectory')
                spe_file = _cache[0] + _cache[2]
                _cache = spe_file.rpartition('WorkingDirectory')
                spe_file = _cache[0] + '/WorkingDirectory/' + _cache[2].split('/')[-1]
            genomes_fullname = id2spe[speid]
            if seq_id == 'nan':
                continue
            record = get_protein(spe_file, seq_id)
            if record is not None:
                record_fullname = id2seq[record.id]
                record.id = genomes_fullname + ' ' + record_fullname
                seqs.append(record)
            else:
                print(seq_id, "doesn't exist?")
        SeqIO.write(seqs, handle=f1, format='fasta-2line')

def get_seq_with_OG(orthogroups_path, OG, output_dir, genomes_list=None,single_copy=True,
                    ):
    os.makedirs(output_dir, exist_ok=True)
    if type(OG) == str:
        OG = [OG]

    thisdir = join(dirname(dirname(orthogroups_path)), 'WorkingDirectory')
    SeqID_file = join(thisdir, "SequenceIDs.txt")
    SpeID_file = join(thisdir, "SpeciesIDs.txt")
    id2seq, seq2id = get_dict(SeqID_file)
    id2spe, spe2id = get_dict(SpeID_file)
    if isinstance(single_copy,pd.DataFrame):
        data = single_copy
    elif isinstance(single_copy,str):
        data = get_single_copy(orthogroups_path)
    else:
        data = pd.read_csv(orthogroups_path, sep='\t', index_col=0)
    if set(OG).difference(set(data.index)):
        raise Exception("Some OG is not presented at the index of data")
    sub_data = data.loc[OG, :]
    if genomes_list is not None:
        #gids = [_ for _ in open(genomes_list).read().split('\n') if _]
        sub_data = data.loc[OG,genomes_list]
    sub_data = sub_data.loc[:, ~sub_data.isna().all(0)]
    # remove all nan genomes
    # open these genomes
    sub_data.columns = [spe2id[_] for _ in sub_data.columns]
    sub_data = sub_data.applymap(lambda x: seq2id.get(str(x).split('.')[0].split(' ')[0], 'nan'))
    species_path_temp = join(thisdir, "Species{speid}.fa")

    params = []
    for og in OG:
        params.append((sub_data,output_dir,og,species_path_temp,id2spe,id2seq))
    with mp.Pool(processes=10) as tp:
        list(tp.imap(extract2file,tqdm(params)))


@click.command()
@click.option("-i", "infile", help='normaly is the file called `Orthogroups.csv`')
@click.option("-og", "OG", help='file or long COMMA separator string,could be None', default=None)
@click.option("-o", "output_dir", help="the directory you want to output to")
@click.option("-single", is_flag=True, default=False, help="You want to use a single copy or not? default is True")
@click.option("-all_g",'contain_all_g',is_flag=True,default=False)
@click.option("-t", "topnumber", help='the minimum number of genomes for each output OG', default=None)
@click.option("-i2", "genomes_list", help="input a file contains the genomes name", default=None)
def main(infile, OG, output_dir, single, topnumber, genomes_list,contain_all_g):
    if OG is None:
        single_copy_df = get_single_copy(infile)
        number_genomes_presence = get_summary_statistic(single_copy_df)
        # give a number to indicate minimum required genomes it need to cover
        if topnumber is not None and genomes_list is not None:
            raise Exception("weird input, -t and -i2 just input one of them. don't specify them both")
        if contain_all_g:
            OG = list(number_genomes_presence.index[number_genomes_presence==single_copy_df.shape[1]])
        if topnumber is not None:
            if OG is not None:
                raise Exception("Overwriting! be careful the input option")
            OG = list(number_genomes_presence.index[:topnumber])

        if genomes_list is not None:
            if OG is not None:
                raise Exception("Overwriting! be careful the input option")
            # give a genomes list file to indicate required genomes
            if not os.path.exists(genomes_list):
                raise Exception("File %s doesn't exist" % genomes_list)
            genomes_list = open(genomes_list, 'r').read().split('\n')
            genomes_list = [_ for _ in genomes_list if _]
            # remove empty string
            weird_names = set(genomes_list).difference(set(single_copy_df.columns))
            if weird_names:
                raise Exception("file you provide contains some weird name, e.g like %s" % '\n'.join(weird_names))
            sub_df = single_copy_df.loc[:, genomes_list]
            OG = list(sub_df.index[~sub_df.isna().any(1)])
    elif os.path.exists(OG):
        OG = open(OG, 'r').read().split('\n')
        OG = [_ for _ in OG if _]
    if OG is None:
        raise Exception("you must specify some reliable parameters, OG now still is NONE")
    get_seq_with_OG(infile, OG, output_dir, genomes_list,single_copy=single)


if __name__ == '__main__':
    main()
    # get single copy of OG which cover at least 200 genomes and output to the current directory
    # python3 ~/script_api/ForOrthofinder/bin/getSeqofOG.py -i /home-user/thliao/project/cyanophage/within_cyanophage/ortho_test/data/Results_Aug09/WorkingDirectory/Orthogroups_3.csv -t 200  -o .
    # python3 ~/script_api/ForOrthofinder/bin/getSeqofOG.py -i 更换成你需要的/Results_Aug09/WorkingDirectory/Orthogroups_3.csv -o 输出目录 -i2 你要指定的基因组的名称的一个list的文件
