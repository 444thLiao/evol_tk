import click
import plotly.graph_objects as go


def coord2polar(x, total, final_degree=360):
    perc = x / total
    theta = final_degree * perc
    return theta

def parse_gff():
    pass

# height with from r0 to dr.
@click.command()
@click.option("-i","infile")

def main():
    pass