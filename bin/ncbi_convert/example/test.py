import os
from subprocess import check_call
from os.path import join

def run(cmd):
    print("testing cmd: ",cmd)
    check_call(cmd,shell=1)
    
current_p = __file__
current_d = os.path.dirname(current_p)
target_d = os.path.dirname(current_d)

data_f = join(current_d,'id_list')

script1 = join(target_d,'pid2GI.py')
script2 = join(target_d,'pid2tax.py')
script3 = join(target_d,'pid2genome.py')

cmd1 = f"python3 {script1} -i {data_f} -o {join(current_d,'test1.tab')} -f"
cmd2 = f"python3 {script2} -i {data_f} -o {join(current_d,'test2.tab')} -f"
cmd3 = f"python3 {script3} -i {data_f} -o {join(current_d,'test3.tab')} -f"

if __name__ == "__main__":
    run(cmd1)
    run(cmd2)
    run(cmd3)