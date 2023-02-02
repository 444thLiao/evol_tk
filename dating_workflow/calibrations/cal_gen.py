import pandas as pd
import itertools
from collections import defaultdict
import sys

def get_mapping(infile='./LCA.tsv'):
    info = open(infile)
    mapping_dict = {}
    for r in info:
        groups,n,l = r.strip().split('\t')
        for g in groups.split(';'):
            mapping_dict[(n,g)] = l
    # from api_tools import pie_chart
    # text = pie_chart({v:{k:1} for k,v in mapping_dict.items()},
    #                     {k:"#ff9900" for k in mapping_dict},)
    # with open('./cal_nodes.txt','w')  as f1:
    #     f1.write(text)    
            
    return mapping_dict


if __name__ == '__main__':
    args = sys.argv[1:]
    
    if len(args)==1:

        cal_info = args[1]
    
    df = pd.read_excel(cal_info, index_col=0)        
    df = df.reindex(columns=mapping_dict)

    setting2node2time = defaultdict(str)
    for s, row in df.iterrows():
        for node, t in row.items():
            if (not pd.isna(t)) and node in mapping_dict :
                setting2node2time[s] += f"{mapping_dict[node]}\t{str(t)}\t{node}\n"
    



    # bac = [f'Bac{i}' for i in '1234']
    # euk = [f'NE{i}' for i in '12345']
    # for c in itertools.product(bac,euk):
    #     final_name = f"{c[0].replace('Bac','B')}{c[1]}"
    #     a = setting2node2time[c[0]]
    #     b = setting2node2time[c[1]]
    #     with open(f"./{final_name}.txt", 'w') as f1:
    #         f1.write(a.strip()+'\n'+b.strip())

# with open(f"./{row.name}.txt", 'w') as f1:
#     f1.write('\n'.join(rows))


