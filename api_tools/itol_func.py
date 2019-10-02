# data dependent transform

def deduced_legend(info2color,info_name='dataset'):
    
    legend_title = info_name
    legend_shape = ','.join(['1'] * len(info2color))
    legend_colors = ','.join([_
                              for _ in info2color.values()])
    legend_labels = ','.join(list(info2color.keys()))
    legend_text = f"""
LEGEND_TITLE,{legend_title}
LEGEND_SHAPES,{legend_shape}
LEGEND_COLORS,{legend_colors}
LEGEND_LABELS,{legend_labels}"""
    return legend_text

def to_color_strip(id2info,info2color,info_name='dataset'):
    template_text = open(
        '/home-user/thliao/template_txt/dataset_color_strip_template.txt').read()
    id2col = {id:info2color[info] for id,info in id2info.items()}
    annotate_text = '\n'.join(['%s,%s\n' % (id,col) for id,col in id2col.items()])
    legend_text = deduced_legend(info2color,info_name)
    
    template_text = template_text.format(legend_text=legend_text,
                     dataset_label=info_name)
    info_name = info_name.replace('/','_')
    return template_text+'\n'+annotate_text
    
def to_color_branch(ID2info,info2color,dataset_name='color branch'):
    # clade for 
    template_text = open(
        '/home-user/thliao/template_txt/dataset_styles_template.txt').read()
    id2col = {ID:info2color[info] for ID,info in ID2info.items()}
    each_template = '{ID}\t{TYPE}\t{WHAT}\t{COLOR}\t{WIDTH_OR_SIZE_FACTOR}\t{STYLE}\t{BACKGROUND_COLOR}\n'
    legend_text = deduced_legend(info2color,dataset_name)
    
    template_text = template_text.format(dataset_label=dataset_name,
                                         legend_text=legend_text)
    
    rows = [each_template.format(ID=ID,
                                 TYPE='branch',
                                 WHAT='node',
                                 COLOR=color,
                                 WIDTH_OR_SIZE_FACTOR=2,
                                 STYLE='normal',
                                 BACKGROUND_COLOR='')
        for ID,color in id2col.items()]
    return template_text + '\n'.join(rows)
        
    