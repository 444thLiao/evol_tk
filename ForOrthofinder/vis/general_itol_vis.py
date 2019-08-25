import pandas as pd
import os, click

base_dir = os.path.dirname(__file__)

family_colors_dict = {"Podoviridae": "#1D2C4D",
                      "Myoviridae": "#2C4066",
                      "Siphoviridae": "#26547C",
                      'Other/Unclassified': "#3594B7"}

host_colors_dict = {"Prochlorococcus": '#71CE37',
                    "Synechococcus": "#FF0000",
                    "unknown": "#7F7C75"}


############################################################
# for color
def get_annotated_text(data_file,
                       id_column,
                       annotate_column,
                       output_file,
                       dataset_label='label1',
                       template_file1=os.path.join(base_dir, 'dataset_color_strip_template.txt'),
                       template_file2=os.path.join(base_dir, 'colors_styles_template.txt'),
                       ):
    template_t = open(template_file1).read()
    template_for_branch_t = open(template_file2).read()
    nuc_df = pd.read_csv(data_file, index_col=None, sep=',')
    if annotate_column == "family":
        l2c = family_colors_dict
    elif annotate_column == "host":
        l2c = host_colors_dict
    else:
        raise Exception
    # labels = list(map(str, set(nuc_df.loc[:, annotate_column])))
    # labels = [_ for _ in labels if _.lower() != 'nan' and _ in l2c]

    legend_title = 'which %s' % annotate_column
    legend_shape = ','.join(['1'] * len(l2c))
    legend_colors = ','.join([l2c[_]
                              for _ in l2c.keys()])
    legend_labels = ','.join(list(l2c.keys()))
    legend_text = """
LEGEND_TITLE,{legend_title}
LEGEND_SHAPES,{legend_shape}
LEGEND_COLORS,{legend_colors}
LEGEND_LABELS,{legend_labels}""".format(legend_title=legend_title,
                                            legend_shape=legend_shape,
                                            legend_colors=legend_colors,
                                            legend_labels=legend_labels)

    annotate_text = ''
    annotate_branch_text = ""
    for rid, row in nuc_df.iterrows():
        if str(row[annotate_column]) not in l2c and annotate_column == 'family':
            id = row[id_column]
            color = l2c['Other/Unclassified']
        elif str(row[annotate_column]) not in l2c and annotate_column == 'host':
            continue
        else:
            id = row[id_column]
            color = l2c[row[annotate_column]]
        annotate_text += ','.join(map(str, [id,
                                            color])) + '\n'
        annotate_branch_text += "{id},clade,{color},normal,2\n".format(id=id,color=color)

    with open(output_file, 'w') as f1:
        f1.write(template_t.format(legend_text=legend_text,
                                   dataset_label=dataset_label) + annotate_text)
    with open(output_file.replace('.txt', '') + '.branch.txt', 'w') as f1:
        f1.write(template_for_branch_t + annotate_branch_text)


############################################################
def relabels_text(data_file,
                  id_column,
                  annotate_column,
                  output_file,
                  template_file=os.path.join(base_dir, 'labels_template.txt'),
                  ):
    nuc_df = pd.read_csv(data_file, index_col=None, sep=',')
    template_t = open(template_file).read()
    annotate_text = ''
    for rid, row in nuc_df.iterrows():
        if str(row[annotate_column]).lower() == "nan":
            continue
        annotate_text += ','.join(map(str, [row[id_column],
                                            row[annotate_column]])) + '\n'
    with open(output_file, 'w') as f1:
        f1.write(template_t + annotate_text)


if __name__ == '__main__':
    t = "C:/Users/thliao/Desktop/cyanophage/annotated_cyanophage.csv"
    get_annotated_text(t,
                       'acc_id',
                       'host',
                       'C:/Users/thliao/Desktop/cyanophage/host_annotated.txt'
                       )
    get_annotated_text(t,
                       'acc_id',
                       'family',
                       'C:/Users/thliao/Desktop/cyanophage/family_annotated.txt'
                       )
    relabels_text(t,
                  'acc_id',
                  'orgname',
                  'C:/Users/thliao/Desktop/cyanophage/labels_reannotated.txt'
                  )
