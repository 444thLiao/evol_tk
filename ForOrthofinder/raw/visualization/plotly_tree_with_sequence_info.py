"""
Script to draw a tree with visualizated blast result
You can showing 155 genome data and its genome difference.
"""

import glob
import os
import re
import pandas as pd
import plotly.graph_objs as go
import plotly
from plotly import tools
# from ..process_gbk.utils import fetch_total_length
### can sepratetly perform extracted_part and whole part blastn result
def construct_plot4_splited(labels_df,fasta_dir,blastx_result,extracted_id,dendrogram_data,draw_loss_nasN = False,draw_only_loss_nasN = False,output_fig = False):
    '''
    big idea: draw each species interval into graph.
    :param labels_df:
    :param fasta_dir:
    :param blastx_result:
    :param element_result:
    :param draw_nasN:
    :param draw_loss_nasN:
    :param draw_only_loss_nasN:
    :param output_fig:
    :return:
    '''
    datas = []
    color = ['rgb(91,155,213)', 'rgb(255,192,0)','rgb(204,51,0)']
    record_sample = []
    ticktext = []
    tickvalues = []
    dict_df = pd.read_csv(labels_df)
    recorded_slider_steps = []
    if type(extracted_id) != list:
        extracted_id = [extracted_id]
    for _i in glob.glob(os.path.join('/home/liaoth/project/NR-Pfam/New_data/extracted_part','*','*.fasta')):
        path = _i
        _i = int(_i.split('/')[-1].split('_')[0])
        #import pdb;pdb.set_trace()
        if 'lossNasAN_gene' in path and os.path.isdir(os.path.dirname(path)+'/backup'):
            datas.append(go.Scatter(x=[_i + 0.5],
                                    y=[-2],
                                    yaxis='y1',
                                    xaxis='x3',
                                    marker=dict(symbol='triangle-down', size=15, color='#F4D00B'),
                                    hoverinfo='text+name',
                                    mode='markers',
                                    name = 'nasN neighbour gene',
                                   ))
        elif 'lossNasAN_gene' in path:
            datas.append(go.Scatter(x=[_i + 0.5],
                                    y=[-2],
                                    yaxis='y1',
                                    xaxis='x3',
                                    marker=dict(symbol='triangle-down', size=15, color='#444'),
                                    hoverinfo='text+name',
                                    mode='markers',
                                    name = 'nasN neighbour gene',
                                   ))

    for _i in extracted_id:
        datas.append(go.Scatter(x=[_i + 0.5],
                                y=[-1],
                                yaxis='y1',
                                xaxis='x3',
                                mode='text',
                                text=str(_i),
                                line=go.Line(color='#000000', width=8),
                                name=_i))

    for path in glob.glob(os.path.join(fasta_dir,'*.fasta')):
        db_id = re.findall('/([0-9]*)_whole_genome', path)[0]
        try:
            db_id_16s = db_id
            db_names = dict_df[dict_df.loc[:,'self_id'] == int(db_id)].loc[:, 'NAME']
        except:
            continue
        contig_num = open(path).read().count('>')
        if contig_num == 1:
            datas.append(
                go.Scatter(x=[0], y=[db_id_16s],
                           mode='markers', marker=dict(symbol='triangle-right', size=15, color='#d62728'),
                           name='complete genome'))
        try:
            input_df = pd.read_csv(os.path.join(blastx_result,'%s_ref_blastx_filtered_65.tab' % str(db_id)), sep='\t')
        except:
            continue
        try:
            contig_name = list(input_df[input_df.blast1 == 4108].index)[0]
            nasN_start = input_df[input_df.blast1 == 4108].iloc[:,5][0]
            nasN_end = input_df[input_df.blast1 == 4108].iloc[:, 6][0]
            nasN_start,nasN_end = min( nasN_start,nasN_end),max( nasN_start,nasN_end)
            if draw_only_loss_nasN:
                continue
        except:
            if not draw_loss_nasN:
                continue
            nasN_start, nasN_end = 0,0
            contig_name = ''
        if contig_name != '':
            record_sample.append(path)

        for _extracted_id in extracted_id:
            opath_format = '/home/liaoth/project/NR-Pfam/New_data/extracted_part_to_whole_blastn/blastn_result/%s_%s_blastn.tab'
            # first is whole genome, second is subject, is extracted part.
            query_path_format = '/home/liaoth/project/NR-Pfam/New_data/extracted_part/{id}/{id}_core_gene.fasta'
            if os.path.isfile(opath_format % (db_id,str(_extracted_id))):
                opath = opath_format % (db_id,str(_extracted_id))
            else:
                continue
            query_path = query_path_format.format(id=_extracted_id)
            if not os.path.isfile(query_path):
                query_path = query_path.replace('_core_gene.fasta','_lossNasAN_gene.fasta')
            total_length = fetch_total_length(query_path)
            try:
                cache = pd.read_csv(opath, sep='\t', header=None)
                cache = cache[cache.loc[:, 3] > 3000]
                if len(cache) == 0:
                    continue
            except:
                continue
            element_contigs = list(cache.index)
            q_starts = list(cache.iloc[:,5])
            q_ends = list(cache.iloc[:,6])
            s_starts = list(cache.iloc[:,7])
            s_ends = list(cache.iloc[:, 8])
            cache_all = list(cache.loc[:, 3])

            for _start,_end in zip(s_starts,s_ends):
                _q_start = q_starts[s_starts.index(_start)]
                _q_end = q_ends[s_starts.index(_start)]
                _q_start, _q_end = min(_q_start, _q_end), max(_q_start, _q_end)
                element_contig = element_contigs[s_starts.index(_start)]
                #print element_contig, contig_name
                if len(set(range(_q_start, _q_end)).intersection(set(range(nasN_start, nasN_end)))) / float(
                                        nasN_end - nasN_start + 1) > 0.95 and element_contig == contig_name and _extracted_id != int(db_id):
                    datas.append(
                        go.Scatter(x=[min(_start,_end)/float(total_length) + _extracted_id,
                                      max(_start,_end)/float(total_length)+ _extracted_id],
                                   y=[int(db_id_16s), int(db_id_16s)],
                                   yaxis='y1',
                                   xaxis='x3',
                                   mode='lines',
                                   line=go.Line(color='#000000', width=8),
                                   name=_extracted_id))

                else:
                    datas.append(go.Scatter(x=[min(_start,_end)/float(total_length)+ _extracted_id,
                                               max(_start,_end)/float(total_length)+ _extracted_id],
                                            y=[int(db_id_16s),int(db_id_16s)],
                                            yaxis='y1',
                                            xaxis='x3',
                                            mode='lines',
                                            line=go.Line(color=color[_extracted_id % len(color)], width=8),
                                            name=_extracted_id,
                                            hoverinfo='none'))
                    recorded_slider_steps.append([len(datas) - 1, cache_all[s_starts.index(_start)]])

        ticktext.append(str(db_names.iloc[0]).partition(' ')[2])
        tickvalues.append(int(db_id_16s))
    record_sample = list(set(record_sample))
    for path in record_sample:
        db_id = re.findall('/([0-9]*)_whole_genome', path)[0]
        input_df = pd.read_csv(os.path.join(blastx_result,'%s_ref_blastx_filtered_65.tab' % str(db_id)),sep='\t')
        input_df = input_df[input_df.blast1 == 4108]
        starts = list(input_df.iloc[:, 7])
        ends = list(input_df.iloc[:, 8])
        total_length = fetch_total_length('/home/liaoth/project/NR-Pfam/protein_db_blast/ref_protein.fasta','4108')
        for _start,_end in zip(starts,ends):
            datas.append(go.Scatter(x=[min(_start,_end)/float(total_length),max(_start,_end)/float(total_length)],
                                    y=[int(db_id),int(db_id)],
                                    yaxis='y1',
                                    xaxis='x2',
                                    mode='lines',
                                    line=go.Line(color='#000000', width=8),
                                    hoverinfo='none'))
    # construct a slider
    # steps = []
    # for length in range(3000,6500,500):
    #     step = dict(method = 'restyle',
    #                 args = ['opacity',[1] * len(datas+dendrogram_data)],
    #                 label = length)
    #     for each_idx,each_length in recorded_slider_steps:
    #         if each_length <= length:
    #             step['args'][1][each_idx] = 0
    #     steps.append(step)
    # sliders = [dict(
    #     active=0,
    #     currentvalue={"prefix": "length: "},
    #     pad={"t": 50},
    #     steps=steps
    # )]

    layout = {'yaxis': {'autorange': 'reversed', 'ticktext': ticktext, 'tickvals': tickvalues},
              #'sliders': sliders,
              'margin': dict(l=250, r=30),
              'showlegend': False,
              'width': 1500 + len(extracted_id)*50,
              'height': 3000,
              }
    fig = tools.make_subplots(rows=1, cols=3, shared_yaxes=True,horizontal_spacing = 0)

    fig.layout.update(layout)
    fig.data.extend(datas)
    fig.data.extend(dendrogram_data)
    fig.layout.yaxis.zeroline = False
    fig.layout.xaxis1.showticklabels = False
    fig.layout.xaxis2.showticklabels = False
    fig.layout.xaxis1.range = [0,0.04]
    fig.layout.xaxis2.range = [0,1]
    fig.layout.xaxis3.range = [min(extracted_id), max(extracted_id)]
    if 1.0/(len(extracted_id)+2) < 0.01:
        fig.layout.xaxis1.domain = [0.0, 0.01]
        fig.layout.xaxis2.domain = [0.01,0.02]
        fig.layout.xaxis3.domain = [0.02, 1.0]
    else:
        fig.layout.xaxis1.domain = [0.0,1.0/(len(extracted_id)+2)]
        fig.layout.xaxis2.domain = [1.0/(len(extracted_id)+2),2.0/(len(extracted_id)+2)]
        fig.layout.xaxis3.domain = [2.0/(len(extracted_id)+2),1.0]
    fig.layout.hovermode = 'closest'
    if not output_fig:
        plotly.offline.plot(fig)
    else:
        return fig