import pandas as pd
def parse_mcmc(infile):
    t = pd.read_csv(infile,sep='\t',index_col=0)
    return t.mean()