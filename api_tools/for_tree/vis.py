from Bio import Phylo

import plotly.graph_objs as go
import plotly
from ete3 import Tree
import io
def get_x_positions(tree):
    """Create a mapping of each clade to its horizontal position.

    Dict of {clade: x-coord}
    """
    depths = tree.depths()
    # If there are no branch lengths, assume unit branch lengths
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    return depths


def get_y_positions(tree):
    """Create a mapping of each clade to its vertical position.

    Dict of {clade: y-coord}.
    Coordinates are negative, and integers for tips.
    """
    maxheight = tree.count_terminals()
    # Rows are defined by the tips
    heights = dict((tip, maxheight - i)
                   for i, tip in enumerate(reversed(tree.get_terminals())))

    # Internal nodes: place at midpoint of children
    def calc_row(clade):
        for subclade in clade:
            if subclade not in heights:
                calc_row(subclade)
        # Closure over heights
        heights[clade] = (heights[clade.clades[0]] +
                          heights[clade.clades[-1]]) / 2.0

    if tree.root.clades:
        calc_row(tree.root)
    return heights


def draw_clade_lines(orientation='horizontal',
                     y_here=0, x_start=0, x_here=0, y_bot=0, y_top=0):
    """Create a line with or without a line collection object.

    Graphical formatting of the lines representing clades in the plot can be
    customized by altering this function.
    """
    if orientation == 'horizontal':
        horizontal_linecollections.append(
            [(x_start, y_here), (x_here, y_here)])
    if orientation == 'vertical':
        vertical_linecollections.append(
            [(x_here, y_bot), (x_here, y_top)])


def draw_clade(clade, x_start, color, lw, fixed_length=None):
    """Recursively draw a tree, down from the given clade."""
    x_here = x_posns[clade]
    y_here = y_posns[clade]
    if clade.is_terminal() and fixed_length is not None:
        x_here = fixed_length
    # Draw a horizontal line from start to here
    draw_clade_lines(orientation='horizontal',
                     y_here=y_here,
                     x_start=x_start,
                     x_here=x_here)

    label = str(clade)
    labels.append(label)
    if clade.clades:
        # for those point non-terminals.
        # Draw a vertical line connecting all children
        y_top = y_posns[clade.clades[0]]
        y_bot = y_posns[clade.clades[-1]]
        # Only apply widths to horizontal lines, like Archaeopteryx
        draw_clade_lines(orientation='vertical',
                         x_here=x_here,
                         y_bot=y_bot,
                         y_top=y_top)

        # Draw descendents
        for child in clade:
            draw_clade(child, x_here, color, lw, fixed_length=fixed_length)


def main(path, fixed_length='max', yscale=1, xscale=1):
    # todo: add curated comments
    global x_posns
    global y_posns
    global horizontal_linecollections
    global vertical_linecollections
    global labels

    if type(path) == Tree:
        tree = Phylo.read(io.StringIO(path.write()), 'newick')
    else:
        tree = Phylo.read(path, 'newick')
    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)
    horizontal_linecollections = []
    vertical_linecollections = []
    labels = []
    if fixed_length == 'max':
        draw_clade(tree.root, 0, 'k', '', fixed_length=2)
        fixed_length = int(sorted([_[0][0] for _ in horizontal_linecollections])[-1]) + 1
        horizontal_linecollections = []
        vertical_linecollections = []
    draw_clade(tree.root, 0, 'k', '', fixed_length=fixed_length)
    datas = []
    labels_draw_text = []
    labels_draw_x = []
    labels_draw_y = []
    labels_backup = labels[::]
    labels = [None
              if 'INT' in _i else _i
              for _i in labels]

    labels[0] = None
    for idx, (_o, _t) in enumerate(horizontal_linecollections):
        horizontal_draws_x = []
        horizontal_draws_y = []
        horizontal_draws_x.append(xscale * _o[0])
        horizontal_draws_y.append(yscale * (_o[1] - 1))
        _x1 = xscale * _t[0]
        _y1 = yscale * (_t[1] - 1)
        if labels[idx]:

            horizontal_draws_x.append(_x1)
            horizontal_draws_y.append(_y1)

            labels_draw_text.append(labels_backup[idx])
            labels_draw_x.append(_x1)
            labels_draw_y.append(_y1)

        else:
            horizontal_draws_x.append(_x1)
            horizontal_draws_y.append(_y1)

            labels_draw_text.append(labels_backup[idx])
            labels_draw_x.append(_x1)
            labels_draw_y.append(_y1)

        trace1 = go.Scatter(x=horizontal_draws_x,
                            y=horizontal_draws_y,
                            mode='lines',
                            line=dict(color='#444', width=1),
                            name=labels_backup[idx],
                            hoverinfo='none', xaxis='x1', yaxis='y1', showlegend=False)
        datas.append(trace1)
    for _o, _t in vertical_linecollections:
        vertical_draws_x = []
        vertical_draws_y = []
        vertical_draws_x.append(xscale * _o[0])
        vertical_draws_x.append(xscale * _t[0])
        vertical_draws_y.append(yscale * (_o[1] - 1))
        vertical_draws_y.append(yscale * (_t[1] - 1))

        trace1 = go.Scatter(x=vertical_draws_x,
                            y=vertical_draws_y,
                            mode='lines',
                            line=dict(color='#444', width=1),
                            hoverinfo='none', xaxis='x1', yaxis='y1', showlegend=False)
        datas.append(trace1)
    return datas, labels, labels_backup, labels_draw_text, labels_draw_x, labels_draw_y


# if __name__ == '__main__':
#     datas, labels, labels_backup, labels_draw_text, labels_draw_x, labels_draw_y = main('./newick.example')
#     layout = dict(yaxis=dict(zeroline=False), width=1400,
#                   height=2000)
#     fig = go.Figure(data=datas, layout=layout)
#     plotly.offline.plot(fig)
