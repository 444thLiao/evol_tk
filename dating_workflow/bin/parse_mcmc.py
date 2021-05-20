import warnings
from glob import glob
from os.path import dirname,exists,join,basename,getsize

import click
import pandas as pd
import plotly.express as px
from ete3 import Tree
from tqdm import tqdm
import re
from dating_workflow.debug_for.draw_infinite_plot import get_plot, fit_line
from dating_workflow.toolkit.mcmctree_for import get_node_name_from_log,get_posterior_df

warnings.filterwarnings("ignore")


def get_node_name(f):
    # f should be the *.out file
    matched_row = ""
    match = False
    for row in open(f):
        row = row.strip("\n").strip()
        if match and not row:
            break
        if row.startswith(
            "Species tree for FigTree.  Branch lengths = posterior mean times; 95% CIs = labels"
        ):
            match = True
            continue
        if match:
            matched_row = row
    t = Tree(matched_row.replace(" ", ""), format=8)
    for l in t.get_leaves():
        l.name = l.name.partition("_")[-1]
    return t


def targetgroup_compare_violin(collect_, odir):
    # todo
    # not finish
    tmp = []
    for k, v in collect_.items():
        _df = pd.DataFrame(v.values)
        _df.columns = ["time"]
        _df.loc[:, "name"] = k
        if "run1" in k or "run" not in k:
            tmp.append(_df)
    df = pd.concat(tmp, axis=0)
    df.loc[:, "num_set"] = [
        int(_.split("_set")[-1].split("_")[0]) if "run" in _ else 0
        for _ in df.loc[:, "name"]
    ]
    df.loc[:, "set"] = [
        _.replace("83g_", "").replace("_run1", "").replace("clock2_", "")
        for _ in df.loc[:, "name"]
    ]
    df = df.sort_values("num_set")

    fig = px.violin(df, x="set", y="time", box=True, points=False)
    fig.layout.yaxis.title.text = "Divergence time(100Mya)"
    fig.layout.xaxis.title.text = "Sets of calibration information"
    fig.layout.yaxis.title.font.size = 30
    fig.layout.xaxis.title.font.size = 30
    fig.layout.xaxis.tickfont.size = 20
    fig.write_html(
        "./dating_for/83g/nucl/set14_derivative.html", include_plotlyjs="cdn"
    )


def main(indir, name2group, odir, no_plot=False, prefix="set", burn_in=0, scale=1):
    tmp_df = pd.DataFrame()
    processed_dir = [
        _ for _ in glob(join(indir, "*run*")) if exists(join(_, "mcmc.txt"))
    ]
    # it would overwrite each when there are run1/run2
    highlight_nodes = []
    for each_dir in tqdm(processed_dir):
        mcmc = join(each_dir, "mcmc.txt")
        if getsize(mcmc) == 0:
            continue
        # outfile = glob(join(each_dir, '*.out'))
        # outfile = [_ for _ in outfile if 'slurm' not in _][0]

        logfile = glob(join(each_dir, "*.log"))
        logfile = [_ for _ in logfile if "slurm" not in _][0]

        set_name = basename(dirname(mcmc)).partition("_")[-1]
        if not set_name:
            set_name = basename(dirname(mcmc))

        t = get_node_name_from_log(logfile)  # get the tree with internal node name
        post_df = get_posterior_df(mcmc, scale=scale, burn_in=burn_in)
        coef, r2 = fit_line(
            x=post_df["Posterior mean time (100 Ma)"].values[:-1],   # remove lnL
            y=post_df["CI_width"].values[:-1]  # remove lnL
        )
        coef = round(coef, 4)
        r2 = abs(round(r2, 4))
        tmp_df.loc[set_name, ["r-square", "slope"]] = r2, coef
        root_name = "t_n%s" % t.name
        tmp_df.loc[set_name, "ROOT"] = "%s (%s) " % (
            post_df.loc[root_name, "Posterior mean time (100 Ma)"],
            post_df.loc[root_name, "CIs"],
        )
        # special
        v = post_df.loc["lnL", "CIs"]
        _t = post_df.loc["lnL", "Posterior mean time (100 Ma)"]
        tmp_df.loc[set_name, "lnL"] = f"{_t} ({v}) "

        for gname, group in name2group.items():
            raw_name = "t_n%s" % t.get_common_ancestor(group).name
            highlight_nodes.append((raw_name, gname))
            tmp_df.loc[set_name, gname] = "%s (%s) " % (
                post_df.loc[raw_name, "Posterior mean time (100 Ma)"],
                post_df.loc[raw_name, "CIs"],
            )
    tmp_df.loc[:, "num_set"] = [
        int(re.findall(f"{prefix}(\d+)", _)[0]) if "run" in _ else 0
        for _ in tmp_df.index
    ]
    final_df = tmp_df.sort_values("num_set")
    # tmp_df.columns = ["divergence time/100Mya (CI)"]

    # draw infinite sites plot
    writer = pd.ExcelWriter(join(odir, f"{basename(odir)}.xlsx"), engine="xlsxwriter")

    pattern = join(indir, "*_run1", "mcmc.txt")
    if no_plot:
        df = get_plot(
            pattern,
            odir,
            no_plot=no_plot,
            highlight_nodes=highlight_nodes,
            prefix=prefix,
        )
        _df = df.copy()
        _df.index = tmp_df.index
        df.to_excel(
            writer, index_label="Calibration sets", sheet_name="infinite site plots"
        )
    final_df.to_excel(
        writer, index_label="Calibration sets", sheet_name="Posterior time"
    )
    writer.save()


@click.command()
@click.option("-i", "indir", help="dir have multiple calibration set")
@click.option("-ns", "target_group", default=None, help="use comma(,) to separate each")
@click.option("-name", "groupname", default="you could separated it with ; ")
@click.option("-disable_plot", "no_plot", default=False, required=False, is_flag=True)
@click.option("-o", "odir", default=None)
@click.option("-p", "prefix", default="set", help="set name prefix, default is set ")
@click.option(
    "-scale",
    "scale",
    default=1,
    help="set scale of the time accordingly, default is 1 ",
)
def cli(indir, target_group, groupname, odir, no_plot, prefix, scale):
    name2group = dict(
        zip(groupname.split(";"), [_.strip() for _ in target_group.split(";")])
    )
    name2group = {k: [_.strip() for _ in v.split(",")] for k, v in name2group.items()}
    # targe_group = [_.strip() for _ in targe_group.split(',')]
    if odir is None:
        odir = indir
    main(
        indir=indir,
        name2group=name2group,
        odir=odir,
        no_plot=no_plot,
        prefix=prefix,
        scale=scale,
    )


if __name__ == "__main__":
    cli()
    # python3 ~/script/evolution_relative/dating_workflow/bin/parse_mcmc.py -i ./dating_for/83g/clock2_diff_cal/ -ns 'GCA_001828545.1,GCA_004282745.1' -name 'Anammox group'

    # ns = ['GCA_001828545.1', 'GCA_004282745.1'] # anammox
    # indir = './dating_for/83g/clock2_diff_cal/'
    # main(indir=indir, ns=ns, groupname='Anammox group', odir=indir)

    # ns = ['GCA_004421105.1', 'GCA_900119085.1'] # AOB
    # indir = './dating/160g/nucl/clock2_diff_cal/'
    # main(indir=indir, ns=ns, groupname='AOB group', odir=indir)
    #
