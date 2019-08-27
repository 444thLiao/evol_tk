import click
import plotly.graph_objects as go


def coord2polar(x,total):
    perc = x/total
    theta = 360 * perc
    return theta

# height with from r0 to dr.

