from api_tools.metadata_for.auto_classify import _classificated,diff_marine_non_marine
from api_tools.IO_for.read import read_table
import click

def main(infile):
    df = read_table(infile,sep='\t',index_col=None)
    pre_df = _classificated(df)
    general_df = diff_marine_non_marine(pre_df)
    general_df.to_excel('./auto_habitat.xlsx',index=0)
    
@click.command()
@click.option("-i")
def cli():
    pass

if __name__ == "__main__":
    cli()