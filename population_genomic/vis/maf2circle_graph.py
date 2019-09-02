import click
from plotly import graph_objs as go
from collections import defaultdict


def coord2polar(x, total, final_degree=360):
    perc = x / total
    theta = final_degree * perc
    return theta


def parse_maf(infile):
    f1 = open(infile, 'r')
    assert next(f1).startswith('##maf')
    usefull = False
    seg2info = defaultdict(dict)
    label = None
    list_db = []
    for line in f1:
        if line.startswith('a'):
            multi = int(line.strip('\n').split('=')[-1])
            label = line.split('=')[-2].split(' ')[0]
            if multi >= 2:
                usefull = True
            else:
                usefull = False
            continue
        if line == '\n' or not usefull:
            continue
        if '\t' in line:
            tag, _src_start, size, strand, srcSIZE, seq = line.strip('\n').split(' ')
            src, start = _src_start.split('\t\t')
            db, contig = (src.partition('.')[0], src.partition('.')[2])
        else:
            tag, src, start, size, strand, srcSIZE, seq = line.strip('\n').split(' ')
            db, contig = (src.partition('.')[0], src.partition('.')[2])
        seg2info[label][db] = (tag, db, contig, start, size, strand, srcSIZE, seq)
        list_db.append(db)
    list_db = sorted(set(list_db))
    return seg2info, list_db


def reorder_with_ref(seg2info, ref=None, total_length=0):
    if all([1 - (ref in v) for v in seg2info.values()]):
        # if all absence
        raise Exception('ref not found in this alignment file.')
    seg_lables = list(seg2info.keys())

    def order_for(label):
        each_infos = seg2info[label]
        if ref not in each_infos:
            return total_length + 1
        else:
            ref_info = each_infos[ref]
            return int(ref_info[3])

    order_labels = sorted(seg_lables, key=order_for)
    return order_labels


def get_pos(seg2info, list_db, ref, polar=True):
    filter_threshold = 2000
    total_length = sum([len(list(segdict.values())[0][-1])
                        for seg, segdict in seg2info.items()
                        if len(list(segdict.values())[0][-1]) >= filter_threshold])
    order_labels = reorder_with_ref(seg2info,
                                    ref=ref,
                                    total_length=total_length)
    plot_data = defaultdict(list)
    r_step = 10
    # only take the length of all segments
    last_end = 0
    for seg_label in order_labels:
        each_infos = seg2info[seg_label]
        info = list(each_infos.values())[0]
        width = len(info[-1])
        if width <= filter_threshold:
            continue
        if polar:
            width_data = coord2polar(width, total_length, final_degree=350)
        else:
            width_data = width/total_length
        half_width_data = width_data / 2
        start_data = half_width_data + last_end
        last_end = half_width_data + start_data
        for _ in list_db:
            # repear multiple times
            plot_data['theta'].append(start_data)
            plot_data['width'].append(half_width_data*2)
        for idx, db in enumerate(list_db):
            if db not in each_infos:
                plot_data['color'].append('white')
            else:
                # todo: add custom color
                plot_data['color'].append('blue')
            plot_data['r'].append(r_step)
    return plot_data


def draw_barpolar(plot_data):
    fig = go.Figure()
    fig.add_trace(go.Barpolar(r=plot_data['r'],
                              width=plot_data['width'],
                              theta=plot_data['theta'],
                              marker=dict(color=plot_data['color'])))
    return fig


def draw_barplot(plot_data):
    fig = go.Figure()
    fig.add_trace(go.Bar(y=plot_data['r'],
                         x=plot_data['width'],
                         marker=dict(color=plot_data['color'])))
    return fig


if __name__ == '__main__':
    pass
#
# # height with from r0 to dr.
# @click.command()
# @click.option("-i", "infile")
# def main():
#     pass
