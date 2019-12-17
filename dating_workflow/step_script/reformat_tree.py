"""
reformat the result tree from formatted into tree without branch and with calibrations.
"""
import click
import os

from api_tools.itol_func import *
from dating_workflow.step_script import process_path

# tree_file = '../trees/iqtree/169g_concat.formatted.newick'
calibration_txt = './calibration.txt'


# root_at = 'GCA_000011385.1,GCA_000009725.1'


def main(tree_file, ofile, itol_odir, calibration_txt=calibration_txt):
    calibration_dict = {}
    for row in open(calibration_txt):
        if row and not row.startswith('#'):
            LCA, time, remained = row.split('\t')
            calibration_dict[LCA] = time

    t = Tree(tree_file, format=3)
    # new_root = t.get_common_ancestor([l 
    #                                   for l in t.get_leaves() 
    #                                   if l.name in root_at])
    # t.set_outgroup(new_root)

    for LCA, time in calibration_dict.items():
        names = LCA.split('|')
        LCA_node = t.get_common_ancestor([l
                                          for l in t.get_leaves()
                                          if l.name in names])
        LCA_node.name = time

    for n in t.traverse():
        if not n.is_leaf():
            if 'I' in n.name and 'S' in n.name:
                n.name = ''
    final_tree = t.write(format=8, format_root_node=True).replace('NoName', '')

    with open(ofile, 'w') as f1:
        f1.write(f'{len(t.get_leaves())}\n' + final_tree)

    size = '12'
    shape = '2'
    filled = '1'
    color = '#ff9900'

    calibration_dict = {}

    rows_str = []
    rows = []
    for row in open(calibration_txt):
        if row and not row.startswith('#'):
            LCA, time, remained = row.split('\t')

            row = '\t'.join([LCA, shape, size, color, filled, '1', time])
            rows.append(row)
            row = '\t'.join([LCA, time, '1', '#FF0000', 'bold', '2', '0'])
            rows_str.append(row)
    template_text = open(dataset_symbol_template).read()
    annotate_text = '\n'.join(rows)
    template_text = template_text.format(dataset_label='calibration',
                                         legend_text='',
                                         maximum_size=size)
    with open(join(itol_odir, 'dating_tree_calibration.txt'), 'w') as f1:
        f1.write(template_text + '\n' + annotate_text)

    template_text = open('/home-user/thliao/template_txt/dataset_text_template.txt').read()
    annotate_text = '\n'.join(rows_str)
    with open(join(itol_odir, 'dating_tree_calibration_str.txt'), 'w') as f1:
        f1.write(template_text + '\n' + annotate_text)

    text = to_node_symbol(tree_file)
    with open(join(itol_odir, 'dating_tree_bootstraps.txt'), 'w') as f1:
        f1.write(text)


@click.command()
@click.option("-i", "intree")
@click.option("-o", "otree")
@click.option("-itol_dir", "itol_dir")
@click.option("-c", "calibration_txt", default=calibration_txt)
def cli(intree, otree, itol_dir, calibration_txt):
    otree = process_path(otree)
    intree = process_path(intree)
    calibration_txt = process_path(calibration_txt)
    itol_dir = process_path(itol_dir)

    if not exists(dirname(otree)):
        os.makedirs(dirname(otree))

    main(tree_file=intree,
         ofile=otree,
         itol_odir=itol_dir,
         calibration_txt=calibration_txt)


if __name__ == "__main__":
    cli()
