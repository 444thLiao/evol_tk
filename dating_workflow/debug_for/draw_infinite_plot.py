import os
from glob import glob
from os.path import *

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from tqdm import tqdm
from dating_workflow.toolkit.mcmctree_for import *

# def get_CI(infile):
#     """

#     :param f: log file?
#     :return:
#     """
#     f = open(infile).read().split('\n')
#     if infile.endswith('.out'):
#         head = 'Posterior mean (95% Equal-tail CI) (95% HPD CI) HPD-CI-width'
#     elif infile.endswith(".log"):
#         head = 'Posterior means (95% Equal-tail CI) (95% HPD CI) HPD-CI-width'
#     if head not in f:
#         return None
#     idx = f.index(head)
#     remained_txt = f[idx + 1:]

#     def format_v(x):
#         x = x.strip('(')
#         x = x.strip(')')
#         x = x.strip(',')
#         return float(x)

#     idx = []
#     CIs = []
#     mean_collect = []
#     CI_collect = []
#     for row in remained_txt:
#         if row.startswith('t_n') or row.startswith('lnL'):
#             vals = row.split(' ')
#             vals = [_ for _ in vals if _ and _ not in '(),']
#             posterior_mean, equal_tailed_95p_down, equal_tailed_95p_up = map(format_v, vals[1:4])
#             CI_collect.append((equal_tailed_95p_up - equal_tailed_95p_down))
#             mean_collect.append(posterior_mean)
#             idx.append(vals[0])
#             CIs.append('%s - %s' % (equal_tailed_95p_down, equal_tailed_95p_up))
#     df = pd.DataFrame()
#     df.loc[:, 'CI_width'] = CI_collect
#     df.loc[:, 'CIs'] = CIs
#     df.loc[:, 'Posterior mean time (100 Ma)'] = mean_collect
#     df.index = idx
#     return df

# def get_r2(text):
#     r_squre_text = text['hovertemplate'].split('<br>')[2]
#     coef = text['hovertemplate'].split('<br>')[1]
#     coef = coef.replace('CI_width', 'y')
#     coef = coef.replace('Posterior mean time (100 Ma)', 'x').split('+')[0].strip()
#     prefix, v = r_squre_text.split('=')

#     return f"{coef}, {prefix} = {v}"


def fit_line(x, y):
    from sklearn.linear_model import LinearRegression

    x = x.reshape(-1, 1)
    y = y.reshape(-1, 1)
    reg = LinearRegression(fit_intercept=False).fit(x, y)

    return reg.coef_[0][0], reg.score(x, y)


def draw_r(df, group=None):
    df = df.copy()
    df = df.reindex([_ for _ in df.index if _.startswith("t_n")])
    df.loc[:, "Posterior mean time (100 Ma)"] = (
        df.loc[:, "Posterior mean time (100 Ma)"] / 10
    )
    df.loc[:, "CI_width"] = df.loc[:, "CI_width"].astype(float) / 10

    fig = px.scatter(
        df, x="Posterior mean time (100 Ma)", y="CI_width", trendline="ols"
    )

    coef, r2 = fit_line(fig.data[0].x, fig.data[0].y)
    coef = round(coef, 4)
    r2 = abs(round(r2, 4))
    # it will generate negative r2, althought it is right ,it may make someone confused
    fig.data[-1].y = coef * fig.data[-1].x

    r_squre_text = ["y = %s * x, R<sup>2</sup> = %s" % (coef, r2)]
    if group is not None:
        xs = list(df.loc[~df[group].isna(), "Posterior mean time (100 Ma)"])
        ys = list(df.loc[~df[group].isna(), "CI_width"])
        ts = list(df.loc[~df[group].isna(), "text"])
        fig.add_scatter(
            x=xs, y=ys, mode="markers", text=ts, marker=dict(color="#ff0000")
        )

    fig.update_layout(
        showlegend=False,
        annotations=[
            go.layout.Annotation(
                x=0.5,
                y=max(list(df.loc[:, "CI_width"])),
                text="<br>".join(r_squre_text),
                showarrow=False,
                font=dict(size=15),
            )
        ],
    )
    fig.layout.xaxis.title.text = "Mean posterior divergence time (Gya)"
    fig.layout.yaxis.title.text = "95% HPD CI width (Gyr)"
    r_squre_v = r2
    return fig, r_squre_v, coef


def get_plot(pattern, odir, no_plot=False, highlight_nodes=None, prefix="set"):
    # pattern should be multiple mcmc.txt
    if not exists(odir):
        os.makedirs(odir)

    tmp_df2 = pd.DataFrame()
    tmp_df = pd.DataFrame()
    a = glob(pattern)
    for f1 in tqdm(a):
        name = f1.split("_")[-3]
        set_name = f1.split("_")[-2]
        df1 = get_posterior_df(f1, burn_in=0, scale=1)
        if df1 is None:
            continue
        if highlight_nodes is not None:
            ns = [_[0] for _ in highlight_nodes]
            ts = [_[1] for _ in highlight_nodes]
            df1.loc[ns, "Highlight"] = "Target"
            df1.loc[ns, "text"] = ts
            # df1 = df1.fillna('')
        fig, r_squre_v, coef = draw_r(df1, group="Highlight")
        if no_plot:
            pass
        else:
            fig.write_image(join(odir, f"repeat_{name}_{set_name}.png"))
        tmp_df.loc[name, set_name] = r_squre_v
        tmp_df2.loc[name, set_name] = coef

    tmp_df = tmp_df.reindex(
        columns=sorted(tmp_df.columns, key=lambda x: int(x.replace(prefix, "")))
    )
    # tmp_df.to_excel(join(odir, 'infinite_site_r2.xlsx'))
    tmp_df2 = tmp_df2.reindex(
        columns=sorted(tmp_df2.columns, key=lambda x: int(x.replace(prefix, "")))
    )
    # tmp_df2.to_excel(join(odir, 'infinite_site_coef.xlsx'))
    new_df = pd.concat([tmp_df.T, tmp_df2.T], axis=1)
    new_df.columns = ["r-square", "slope"]
    new_df.to_excel(join(odir, "infinite_site.xlsx"), index_label="Calibration sets")
    return new_df


def separate_fit(df, odir, id_sets):
    total_df = df.copy()
    for idx, _ in enumerate(id_sets):
        total_df.loc[_, "color"] = [f"set{idx + 1}"]
    total_df = total_df.fillna("set")
    fig, r_squre_v, coef = draw_r(total_df, group="color")

    _df_sets = [df.loc[_, :].copy() for _ in id_sets]

    results = [draw_r(_df) for _df in _df_sets]

    for idx, (_f, _coef, _r2) in enumerate(results):
        _data = [
            _
            for _ in fig.data
            if str(idx + 1) in _["legendgroup"]
            and "OLS trendline" in _["hovertemplate"]
        ][0]
        _data.y = _f.data[-1].y
        _data.x = _f.data[-1].x

    df.loc[:, "Posterior mean time (100 Ma)"] = (
        df.loc[:, "Posterior mean time (100 Ma)"] / 10
    )
    df.loc[:, "CI_width"] = df.loc[:, "CI_width"].astype(float) / 10
    fig = px.scatter(
        df, x="Posterior mean time (100 Ma)", y="CI_width", trendline="ols"
    )
    coef, r2 = fit_line(fig.data[0].x, fig.data[0].y)
    coef = round(coef, 4)
    r2 = abs(round(r2, 4))
    # it will generate negative r2, althought it is right ,it may make someone confused
    fig.data[2].y = coef * fig.data[2].x

    r_squre_text = ["y = %s * x, R<sup>2</sup> = %s" % (coef, r2)]

    fig.update_layout(
        showlegend=False,
        annotations=[
            go.layout.Annotation(
                x=0.5,
                y=max(list(df.loc[:, "CI_width"])),
                text="<br>".join(r_squre_text),
                showarrow=False,
                font=dict(size=15),
            )
        ],
    )
    fig.layout.xaxis.title.text = "Mean posterior divergence time (Gya)"
    fig.layout.yaxis.title.text = "95% HPD CI width (Gyr)"
    r_squre_v = r2
    return fig, r_squre_v, coef


if __name__ == "__main__":
    odir = "./dating_for/83g/clock2_infinite_plot"
    pattern = "./dating_for/83g/clock2_diff_cal/*_run1/run.log"
    get_plot(pattern, odir)

    plancto_set = [f"t_n{num}" for num in range(139, 166)]
    cyano_set = [f"t_n{num}" for num in range(85, 139)]
    id_sets = [plancto_set, cyano_set]
    odir = "./dating_for/83g/clock2_infinite_plot"
    pattern = "./dating_for/83g/clock2_diff_cal/*_run1/run.log"
    get_plot(pattern, odir)
