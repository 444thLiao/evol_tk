import click
from plotly import graph_objs as go
from collections import defaultdict
import pandas as pd


def coord2polar(x, total, final_degree=360):
    perc = x / total
    theta = final_degree * perc
    return theta


def parse_OG_tab(infile, subset_samples=None,final_degree=350):
    OG_df = pd.read_csv(infile, sep='\t', index_col=0)
    if subset_samples is not None:
        subset_samples = open(subset_samples).read().split('\n')
        OG_df = OG_df.loc[:,subset_samples]
        OG_df = OG_df.loc[~OG_df.isna().all(1),:]
    num_rows = OG_df.shape[0]
    plot_data = defaultdict(list)
    width = final_degree / num_rows
    half_width = width / 2
    last_end = 0
    r_step = 10
    for OG_id, row in OG_df.iterrows():
        _colors = []
        theta = half_width + last_end
        last_end = half_width + theta
        for locus in row.values:
            plot_data['theta'].append(theta)
            plot_data['width'].append(width)
            if pd.isna(locus):
                plot_data['color'].append('white')
            else:
                plot_data['color'].append('red')
            plot_data['r'].append(r_step)
    return plot_data


def draw_barpolar(plot_data):
    fig = go.Figure()
    fig.add_trace(go.Barpolar(r=plot_data['r'],
                              width=plot_data['width'],
                              theta=plot_data['theta'],
                              marker=dict(color=plot_data['color'])))
    fig.update_layout(
        showlegend=False,
        polar=dict(
            bgcolor="white",
            angularaxis=dict(
                visible=False, ),
            radialaxis=dict(
                visible=False)
        ),
        paper_bgcolor="white",
        plot_bgcolor='white'

    )
    return fig


def draw_barplot(plot_data):
    fig = go.Figure()
    fig.add_trace(go.Bar(y=plot_data['r'],
                         x=plot_data['width'],
                         marker=dict(color=plot_data['color'])))
    return fig


@click.command()
@click.option("-i", "input_OG")
@click.option("-o", "output_file")
@click.option("-s", "subset_samples")
def cli(input_OG, output_file, subset_samples):
    plot_data = parse_OG_tab(input_OG)


if __name__ == '__main__':
    cli()
