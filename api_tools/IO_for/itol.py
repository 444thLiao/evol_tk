
from itolapi import Itol
import os
apikey = os.environ.get('ITOLKEY')
def get_itoltree(tpath,outfile='./tmp.png',name='tmp',anno_files=[]):
    os.system(f"cp -r {tpath} ./tmp.tree ")
    itol_uploader = Itol()
    itol_uploader.params['projectName'] = 'batch_access_th'
    itol_uploader.params['APIkey'] = apikey
    itol_uploader.params['treeName'] = name
    itol_uploader.add_file('./tmp.tree')
    for f in anno_files:
        itol_uploader.add_file(f)
    status = itol_uploader.upload()
    assert status != False
    itol_exporter = itol_uploader.get_itol_export()
    for _ in range(len(anno_files)):
        itol_exporter.set_export_param_value('datasets_visible',str(_))
    itol_exporter.set_export_param_value('display_mode','2')
    itol_exporter.set_export_param_value('range_mode','2')
    itol_exporter.set_export_param_value('dashed_lines','0')
    suffix = outfile.split('.')[-1]
    itol_exporter.set_export_param_value('format', suffix)
    itol_exporter.export(outfile)
    os.system(f"rm ./tmp.tree")
    return itol_uploader.get_webpage()