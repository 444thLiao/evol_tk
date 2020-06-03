import plotly.express as px
from tqdm import tqdm
import pandas as pd

# assert plotly.__version__ == '4.8.1'


def rename_ambiguous(complete_df):
    # process the copy one instead of the original one
    complete_df = complete_df.copy()
    level2names = {l: complete_df[l].unique() for l in complete_df.columns}  # besides species
    dup_names = []
    for l, names in level2names.items():
        other_l = set(level2names).difference({l})
        vals = [n for ol in other_l for n in level2names[ol]]
        for n in names:
            if pd.isna(n):
                continue
            if n in vals:
                dup_names.append(n)
    dup_names = set(dup_names)
    if dup_names:
        print(f"detect {len(dup_names)} ambiguous names")
    for dn in dup_names:
        for col in complete_df.columns:
            if dn in complete_df[col].values:
                complete_df.loc[complete_df[col] == dn, col] = f"{dn}({col[0]})"
    return complete_df


def generate_sunburst(df):
    """
    df should look like
                          phylum                  class              order              family           genus
GCA_006184665     Proteobacteria    Gammaproteobacteria   Enterobacterales  Enterobacteriaceae      Salmonella
GCA_002384685     Proteobacteria     Betaproteobacteria   Nitrosomonadales                 NaN             NaN
GCA_008284035     Proteobacteria  Epsilonproteobacteria  Campylobacterales  Campylobacteraceae   Campylobacter
GCA_008965025     Proteobacteria    Gammaproteobacteria   Enterobacterales  Enterobacteriaceae      Salmonella
GCA_005192855     Proteobacteria  Epsilonproteobacteria  Campylobacterales  Campylobacteraceae   Campylobacter
    The name of index isn't matter. The column could be other things...
    But the order of the columns would decide the hierarchical relation of the final graph.
    :param df:
    :return:
    """
    data = dict(
        character=[],
        parent=[],
        value=[],
    )
    for idx, col in tqdm(enumerate(df.columns)):
        if idx != 0:
            pre_col = df.columns[idx - 1]
            uniq_vals = [(_1, _2) for _1, _2 in df.loc[:, [pre_col, col]].drop_duplicates().values
                         if _2]
        else:
            _uv = df[col].unique()
            uniq_vals = zip([''] * len(_uv), _uv)

        for pre_v, v in uniq_vals:
            val = df.loc[df[col] == v, :].shape[0]
            data['character'].append(v)
            data['parent'].append(pre_v)
            data["value"].append(val)

    fig = px.sunburst(
        data,
        names='character',
        parents='parent',
        values='value',
        branchvalues="total",
    )
    return fig
