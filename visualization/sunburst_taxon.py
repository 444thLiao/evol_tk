import pandas as pd
import plotly.express as px
from tqdm import tqdm


# assert plotly.__version__ == '4.8.1'
def filled_unassigned(df):
    """
    fill unassigned level using parental name with U. prefixed
    """
    df = df.copy()
    df.phylum = df.phylum.fillna("Unknown Phylum")
    for idx, level in enumerate(df.columns[1:-1]):
        sub_df = df.loc[df.loc[:, level].isna(), :]
        previous_names = ["U. " + str(_) for _ in sub_df.loc[:, df.columns[idx]]]
        df.loc[df.loc[:, level].isna(), level] = previous_names
    return df

#
# def check_ambiguous_pre(complete_df):
#     complete_df = complete_df.copy()
#     level2names = {l: complete_df[l].unique()
#                    for l in complete_df.columns}  # besides species
#
#     levele2ambiguous_names_info = {}
#     for l, names in tqdm(level2names.items()):
#         if l == 'phylum':
#             continue
#         pre_l = list(complete_df.columns).index(l) - 1
#         pre_l = list(complete_df.columns)[pre_l]
#         for name in tqdm(names):
#             u = complete_df.loc[complete_df[l] == name, pre_l].unique()
#             if len(u) != 1:
#                 # print(name,u)
#                 levele2ambiguous_names_info[l] = (name, pre_l, u)


def rename_ambiguous(complete_df):
    # process the copy one instead of the original one
    complete_df = complete_df.copy()
    level2names = {l: complete_df[l].unique()
                   for l in complete_df.columns}  # besides species

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
            uniq_vals = [(_1, _2)
                         for _1, _2 in df.loc[:, [pre_col, col]].drop_duplicates().values
                         if _2]
            # here I filter out the no defined names in this level
        else:
            pre_col = ''
            _uv = df[col].unique()
            uniq_vals = zip([''] * len(_uv), _uv)

        for pre_v, v in uniq_vals:
            if idx != 0:
                val = df.loc[(df[col] == v) & (df[pre_col] == pre_v), :].shape[0]
            else:
                val = df.loc[(df[col] == v), :].shape[0]
            
            if idx != len(df.columns)-1:
                v = v + f"count: {val}"
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


if __name__ == '__main__':
    tax_tab2 = "/home-user/thliao/.cache/ncbi-genome-download/bacteria2taxonomy.tab"
    sub_df2 = pd.read_csv(tax_tab2, sep='\t', index_col=0)
    sub_df2_renamed = rename_ambiguous(sub_df2)
    sub_df2_renamed_filled = filled_unassigned(sub_df2_renamed)
    fig = generate_sunburst(sub_df2_renamed_filled.iloc[:, :-1])

    #
    # count = 0
    # c2v = {}
    # for c,p,v in zip(fig.data[0]['labels'],
    #                     fig.data[0]['parents'],
    #                  fig.data[0]['values']):
    #     if p == 'Candidatus Babeliae':
    #         # print(c,v)
    #         count +=v
    #         c2v[c]=v
    #
    # c2v2 = sub_df2_renamed_filled.loc[sub_df2_renamed_filled['class']=='Candidatus Babeliae'].groupby('order').size()
    # for c in c2v:
    #     if c2v[c] != c2v2[c]:
    #         print(c,c2v[c] , c2v2[c])
