import os,sys
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import click
from itol_func import *
from os.path import *
import pandas as pd
def parse_table(input_table):
    try:
        in_df = pd.read_csv(input_table,sep='\t',index_col=0)
    except:
        in_df = pd.read_csv(input_table,sep=',',index_col=0)
    # todo, define a standard for input table
        

@click.command()
@click.option("-i","input_table",help='a table which contains infos')
@click.option("-t","tree_file",help='tree file',required=False,default=None)
@click.option("-O","output_dir")
@click.option("-m","mode",
              type=click.Choice(['clade_branch', 'cb'],
                                 case_sensitive=False),
              help='cb mean `clade_branch` which colorized branch within whole clade with define standard.')
def cli(input_table,tree_file,output_dir,mode):
    # assert input/output directory
    output_dir = os.path.realpath(os.path.abspath(output_dir))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(input_table):
        raise IOError(f'{input_table} does not exists')
    # if tree_file, renamed it.
    if tree_file is not None:
        out_treefile = join(output_dir,basename(tree_file))
        new_tree = renamed_tree(tree_file,out_treefile)
    # if mode is clade branch
    

if __name__ == "__main__":
    cli()
