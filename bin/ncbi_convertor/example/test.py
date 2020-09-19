import os
from os.path import join
from subprocess import check_call
from bin.ncbi_convertor import NCBI_convertor

def run(cmd):
    print("testing cmd: ", cmd)
    check_call(cmd, shell=1)


current_p = __file__
current_d = os.path.dirname(current_p)
target_d = os.path.dirname(current_d)

data_f = join(current_d, 'protein_ids')

script1 = join(target_d, 'pid2GI.py')
script2 = join(target_d, 'pid2tax.py')
script3 = join(target_d, 'pid2genome.py')
script4 = join(target_d, 'pid2bio.py')

cmd1 = f"python3 {script1} -i {data_f} -o {join(current_d, 'test1.tab')} -f"
cmd2 = f"python3 {script2} -i {data_f} -o {join(current_d, 'test2.tab')} -f"
cmd3 = f"python3 {script3} -i {data_f} -o {join(current_d, 'test3.tab')} -f"
cmd4 = f"python3 {script4} -i {data_f} -o {join(current_d, 'test4.tab')} -f"



# test
pids = open('./protein_ids').read().split('\n')
convertor = NCBI_convertor(pids, db='protein')
# convertor.check_cache(suffix=suffix, redo=redo)
convertor.get_taxons_from_tid()
pid2assembly_dict = convertor.get_protein_pos_INFO()
aids = open('./assembly_ids').read().split('\n')
convertor = NCBI_convertor(aids, db='assembly')
# convertor.check_cache(suffix=suffix, redo=redo)
convertor.get_taxons_from_tid()


if __name__ == "__main__":
    run(cmd1)
    run(cmd2)
    run(cmd3)
    run(cmd4)
