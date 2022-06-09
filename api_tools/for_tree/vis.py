from Bio import Phylo

import plotly.graph_objs as go
from ete3 import Tree
import io


class tree_vis(object):
    def __init__(self, newick, branch_length=True, leaves2top=True):
        # self.newick = newick
        self.tree_obj = self.read_newick(newick)
        self.max_depth = self.get_x_positions(get_max=True, branch_length=branch_length)
        self.leaves2top = leaves2top
        self.vertical_linecollections = []
        self.horizontal_linecollections = []
        # pair of coordinates form a line.
        self.clade_x_pos = self.get_x_positions(
            branch_length=branch_length, leaves2top=leaves2top
        )
        self.clade_y_pos = self.get_y_positions()
        self.labels = []
        self.labels_info = []
        self.idx2clade = {}
        self.yscaled = 1
        self.xscaled = 1
    def refresh(self, branch_length=True, leaves2top=True):
        self.clade_x_pos = self.get_x_positions(
            branch_length=branch_length, leaves2top=leaves2top
        )

    def read_newick(self, newick):
        if type(newick) == Tree:
            return Phylo.read(io.StringIO(newick.write()), "newick")
        elif 'tree_obj' in dir(newick):
            return newick.tree_obj
        else:
            return Phylo.read(newick, "newick")

    @property
    def root(self):
        return self.tree_obj.root

    def get_fig(self,interval=1.5,width=600,height=1200):
        # generate a suitable figure which has automodify the layout including the template, size, ticks
        fig = go.Figure()
        fig.update_layout(
            xaxis=dict(
                zeroline=False,
              
            ),
            yaxis=dict(showticklabels=False, zerolinecolor="#000000"),
            plot_bgcolor="#ffffff",
            # template='simple_white',
            width=width,
            height=height,
        )
        if self.leaves2top:
            v = int(self.max_depth//interval)
            fig.update_layout(     
                              xaxis=dict(          
                              tickmode = 'array',
                              
                tickvals = [-self.max_depth + (interval*i) for i in range(v) ],
                ticktext = [round(self.max_depth- (interval*i),2) for i in range(v)]  ))
        return fig

    def get_index(self,l1,l2,whole_clade=True):
        # get the indexs of lines with particular name
        
        # if you want to annotate particular clade, you could retrieve the index of datas which are the lines descending to the clade
  
        clade = self.tree_obj.common_ancestor([l1,l2])
        idx_collections = []
        for idx,_clade in self.idx2clade.items():
            if clade.is_parent_of(_clade):
                idx_collections.append(idx)
        return idx_collections
    
    def clade2name(self,clade):
        return clade.get_terminals()[0].name + '|' + clade.get_terminals()[-1].name
    
    def get_clade(self,label):
        if '|' in label:
            l1,l2 = label.split('|')
            return self.tree_obj.common_ancestor([l1,l2])
        else:
            return [_ for _ in self.tree_obj.get_terminals() if _.name==label][0]
        
    def get_plotly_data(self, xscale=1, yscale=1, color="#000000", width=1, x_shift=0, y_shift=0, fix_length=None,
    ):
        """[summary]

        Args:
            xscale (int, optional): [description]. Defaults to 1.
            yscale (int, optional): [description]. Defaults to 1.
            color ([type], optional): [description]. Defaults to None.
            x_shift (int, optional): positive mean shift it to right. vice versus. It is applied before the scale. Defaults to 0.
            y_shift (int, optional): positive mean shift it to above. vice versus. It is applied before the scale. Defaults to 0.

        Returns:
            [type]: [description]
        """
        if self.leaves2top:
            root_pos = 0 - self.max_depth
            self.draw_clade(self.root, root_pos, "k", "",fix_length)
        else:
            self.draw_clade(self.root, 0, "k", "",fix_length)
        self.yscaled = yscale
        self.xscaled = xscale
        
        datas = []
        labels = [None if "INT" in _i else _i for _i in self.labels]
        labels[0] = None
        
        datas_idx = len(datas)
        for idx, (left_coord, right_coord,clade) in enumerate(
            self.horizontal_linecollections
        ):
            horizontal_draws_x = []
            horizontal_draws_y = []

            scaled_left_x = xscale * (left_coord[0] + x_shift)
            scaled_left_y = yscale * (left_coord[1]  + y_shift)
            scaled_right_x = xscale * (right_coord[0] + x_shift)
            scaled_right_y = yscale * (right_coord[1]  + y_shift)
            horizontal_draws_x.extend([scaled_left_x, scaled_right_x])
            horizontal_draws_y.extend([scaled_left_y, scaled_right_y])
            if labels[idx] is not None:
                self.labels_info.append([scaled_left_x,
                                         scaled_left_y,
                                         self.labels[idx]]) 

            trace1 = go.Scatter(
                x=horizontal_draws_x,
                y=horizontal_draws_y,
                mode="lines",
                line=dict(color=color, width=width),
                name=self.labels[idx],
                hoverinfo="none",
                xaxis="x1",
                yaxis="y1",
                showlegend=False,
            )
            datas.append(trace1)
            self.idx2clade[datas_idx] = clade
            datas_idx+=1
        for down_coord, above_coord,clade in self.vertical_linecollections:
            vertical_draws_x = []
            vertical_draws_y = []
            scaled_down_x = xscale * (down_coord[0] + x_shift)
            scaled_down_y = yscale * (down_coord[1]  + y_shift)
            scaled_above_x = xscale * (above_coord[0] + x_shift)
            scaled_above_y = yscale * (above_coord[1]  + y_shift)
            vertical_draws_x.extend([scaled_above_x, scaled_down_x])
            vertical_draws_y.extend([scaled_above_y, scaled_down_y])

            trace1 = go.Scatter(
                x=vertical_draws_x,
                y=vertical_draws_y,
                mode="lines",
                line=dict(color=color, width=width),
                hoverinfo="none",
                xaxis="x1",
                yaxis="y1",
                showlegend=False,
            )
            datas.append(trace1)
            self.idx2clade[datas_idx] = clade
            datas_idx+=1            
        return datas

    def add_marker(self,l1,l2,style={'size':10}):
        clade = self.tree_obj.common_ancestor([l1,l2])
        y_pos = self.clade_y_pos[clade]        
        x_here = self.clade_x_pos[clade]

        data = go.Scatter(x=[x_here*self.xscaled],y=[(y_pos)*self.yscaled],
                          mode='markers',marker=style)
        return data
    
    def draw_clade(self, clade, x_start, color, lw, fixed_length=None):
        """Recursively draw a tree, down from the given clade."""
        x_here = self.clade_x_pos[clade]
        y_here = self.clade_y_pos[clade]

        if clade.is_terminal() and fixed_length is not None:
            x_here = fixed_length

        # Draw a horizontal line from start to here
        self.add_line(
            orientation="horizontal", y_here=y_here, x_start=x_start, x_here=x_here,clade=clade
        )

        label = str(clade)
        if label== 'Clade':
            label = clade.get_terminals()[0].name + '|' + clade.get_terminals()[-1].name
        self.labels.append(label)

        if clade.clades:
            # for those point which are non-terminals.
            # Draw a vertical line connecting all children
            y_top = self.clade_y_pos[clade.clades[0]]
            y_bot = self.clade_y_pos[clade.clades[-1]]
            # Only apply widths to horizontal lines
            self.add_line(
                orientation="vertical", x_here=x_here, y_bot=y_bot, y_top=y_top,clade=clade
            )

            # Draw descendents
            for child in clade:
                self.draw_clade(child, x_here, color, lw, fixed_length=fixed_length)

    def add_line(
        self, orientation="horizontal", y_here=0, x_start=0, x_here=0, y_bot=0, y_top=0,clade=None
    ):
        """Create a line with or without a line collection object.

        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if orientation == "horizontal":
            self.horizontal_linecollections.append(
                [(x_start, y_here), (x_here, y_here),clade]
            )
        if orientation == "vertical":
            self.vertical_linecollections.append([(x_here, y_bot), (x_here, y_top),clade])

    def get_y_positions(self):
        """Create a mapping of each clade to its vertical position.

        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        maxheight = self.tree_obj.count_terminals()
        # Rows are defined by the tips
        heights = dict(
            (tip, maxheight - i)
            for i, tip in enumerate(reversed(self.tree_obj.get_terminals()))
        )

        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (
                heights[clade.clades[0]] + heights[clade.clades[-1]]
            ) / 2.0

        if self.tree_obj.root.clades:
            calc_row(self.tree_obj.root)
        return heights

    def get_x_positions(self, branch_length=False, leaves2top=False, get_max=False):
        """Create a mapping of each clade to its horizontal position.

        Dict of {clade: x-coord}
        """
        if branch_length:
            x_posns = self.tree_obj.depths()
        else:
            x_posns = self.tree_obj.depths(unit_branch_lengths=True)
        max_d = max(x_posns.values())
        if get_max:
            # return the total depth of the tree(positive value)
            return max_d
        if leaves2top:
            # put the leaves on the 0, the orther would be negative
            x_posns = {c: d - max_d for c, d in x_posns.items()}
        else:
            # put the root on the 0
            pass

        return x_posns


