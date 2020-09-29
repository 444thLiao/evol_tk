"""
Primarily dependent on the visualization of kepler-mapper
# todo:...

"""

import networkx as nx
import kmapper as km
import os
from jinja2 import Environment, FileSystemLoader, Template

module_root = os.path.join(os.path.dirname(__file__), "templates")
env = Environment(loader=FileSystemLoader(module_root))

def vis():
    color_function = init_color_function(graph, None)

    X_names = []
    lens_names = []
    lens = None
    custom_tooltips = None
    nbins = 10
    title = "Kepler Mapper"
    X = None

    mapper_data = format_mapper_data(
        graph,
        color_function,
        X,
        X_names,
        lens,
        lens_names,
        custom_tooltips,
        env,
        nbins,
    )

    colorscale = colorscale_default

    histogram = graph_data_distribution(graph, color_function, colorscale)

    mapper_summary = format_meta(graph, custom_meta)

    # Find the absolute module path and the static files
    js_path = os.path.join(os.path.dirname(__file__), "static", "kmapper.js")
    with open(js_path, "r") as f:
        js_text = f.read()

    css_path = os.path.join(os.path.dirname(__file__), "static", "style.css")
    with open(css_path, "r") as f:
        css_text = f.read()

    # Render the Jinja template, filling fields as appropriate
    template = env.get_template("base.html").render(
        title=title,
        mapper_summary=mapper_summary,
        histogram=histogram,
        dist_label="Node",
        mapper_data=mapper_data,
        colorscale=colorscale,
        js_text=js_text,
        css_text=css_text,
        show_tooltips=True,
    )

    if save_file:
        with open(path_html, "wb") as outfile:
            if self.verbose > 0:
                print("Wrote visualization to: %s" % (path_html))
            outfile.write(template.encode("utf-8"))



def nxG2mapperG():
    pass