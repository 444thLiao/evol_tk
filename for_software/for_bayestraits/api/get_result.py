import pandas as pd
import io
from collections import defaultdict
infile = 'm2nm.txt.Log.txt'

rows = open(infile).readlines()
header = [(idx,_) for idx,_ in enumerate(rows) if _.startswith('Iteration')][1]
start_at = header[0]
headers = header[1].strip('\n')
headers = headers.split('\t')

result_df = pd.read_csv(io.StringIO(''.join(rows[start_at:]) ),sep='\t')

mean_vals = result_df.mean()

mean_v = mean_vals.to_dict()

transition_rate = {k:v
                   for k,v in mean_v.items() if k.startswith('q')}

n2cat2prob = defaultdict(lambda :defaultdict(dict))
for key,v in mean_v.items():
    if ' P(' in key:
        node_name = key.split(' ')[0]
        if node_name == 'Root':
            node_name = 'OROOT' # for itol
        cat = key.split('(')[-1].strip(')')
        n2cat2prob[node_name][cat] = v
        
cat2info = {"M":'#0000ff',
            "N":'#D68529'}
from api_tools.itol_func import *
text = pie_chart(n2cat2prob,cat2info,)