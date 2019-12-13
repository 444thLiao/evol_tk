import pandas as pd
import io
infile = 'm2nm.txt.Log.txt'

rows = open(infile).readlines()
header = [(idx,_) for idx,_ in enumerate(rows) if _.startswith('Iteration')][1]
start_at = header[0]
headers = header[1].strip('\n')
headers = headers.split('\t')

result_df = pd.read_csv(io.StringIO(''.join(rows[start_at:]) ),sep='\t')

mean_vals = result_df.mean()