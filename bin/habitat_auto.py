from api_tools.metadata_for.auto_classify import *
from api_tools.IO_for.read import read_table
import click

def main(infile):
    df = read_table(infile,sep='\t',index_col=None)
    pre_df = _classificated(df)

@click.command()
@click.option("-i")
def cli():
    pass

if __name__ == "__main__":
    cli()