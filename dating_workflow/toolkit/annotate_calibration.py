from api_tools.itol_func import *

size = '12'
shape = '2'
filled = '1'
color = '#ff9900'

calibration_txt = './calibration.txt'
calibration_dict = {}
rows = []
for row in open(calibration_txt):
    if row and not row.startswith('#'):
        LCA,time,remained = row.split('\t')
        
        row = '\t'.join([LCA,shape,size,color,filled,'1',time])
        rows.append(row)
        
template_text = open(dataset_symbol_template).read()
annotate_text = '\n'.join(rows)
template_text = template_text.format(dataset_label='calibration',
                                        legend_text='',
                                        maximum_size=size)
with open('./itol_txt/itol_calibration.txt','w') as f1:
    f1.write(template_text+'\n'+annotate_text)