import pandas as pd
from collections import Counter
_df = pd.read_excel('./1.xlsx',header=None,index_col=None)
a = list(_df.iloc[:,0])
b = list(_df.iloc[:,0])
df = pd.read_csv(io.StringIO('\n'.join(['\t'.join(_) for _ in a])),sep='\t',header=None)
for _,col in df.iteritems():
    print(sorted({k:round(v/df.shape[0]*100,2) for k,v in Counter(col).items()}.items()))
df = pd.read_csv(io.StringIO('\n'.join(['\t'.join(_) for _ in b])),sep='\t',header=None)
for _,col in df.iteritems():
    print(sorted({k:round(v/df.shape[0]*100,2) for k,v in Counter(col).items()}.items()))