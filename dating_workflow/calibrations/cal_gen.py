import pandas as pd
import itertools
from collections import defaultdict
import sys

def get_mapping(infile='./LCA.tsv'):
    info = open(infile)
    mapping_dict = {}
    for r in info:
        if r.startswith('#'):continue
        groups,n,l = r.strip().split('\t')
        for g in groups.split(';'):
            mapping_dict[(n.split(' ')[0],g)] = (l,n)
    return mapping_dict

if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args)==1:
        cal_info = args[1]
    else:
        exit('no cal_info found')
    mapping_dict = get_mapping()
    
    from api_tools import pie_chart
    text = pie_chart({v[0]:{v[1]:1} for k,v in mapping_dict.items()},
                     {v[1]:"#ff9900" for k,v in mapping_dict.items()},)
    with open('./cal_nodes.txt','w')  as f1:
        f1.write(text)   
    from api_tools.itol_func import dataset_text_template
    
    template_text = open(dataset_text_template).read()
    annotate_text = '\n'.join(rows_str)
    with open(join(itol_odir, 'dating_tree_calibration_str.txt'), 'w') as f1:
        f1.write(template_text + '\n' + annotate_text)
         
    df = pd.read_excel(cal_info, index_col=0)     
    df.columns = [_.split(' ')[0] for _ in df.columns]   
    df = df.reindex(columns=[_ for _ in mapping_dict])

    setting2node2time = defaultdict(str)
    for s, row in df.iterrows():
        for node, t in row.items():
            if (not pd.isna(t)) and node in mapping_dict:
                info = mapping_dict[node]
                setting2node2time[s] += f"{info[0]}\t{str(t)}\t{info[1]}\n"
    # all_sets = ['B','A','P']
    # set2cals = {i:[_ for _ in setting2node2time if _.startswith(i)] for i in all_sets}
    # for i in set2cals['B']:
    #     with open(f"/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/system_dating/cal/{i}.txt", 'w') as f1:
    #         f1.write(setting2node2time[i])        
    # for b1,a1 in itertools.product(set2cals['B'],set2cals['A']):
    #     with open(f"/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/system_dating/cal/{b1}{a1}.txt", 'w') as f1:
    #         f1.write(setting2node2time[b1]+'\n'+setting2node2time[a1])  
    # for b1,a1 in itertools.product(set2cals['B'],set2cals['P']):
    #     with open(f"/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/system_dating/cal/{b1}{a1}.txt", 'w') as f1:
    #         f1.write(setting2node2time[b1]+'\n'+setting2node2time[a1])  


