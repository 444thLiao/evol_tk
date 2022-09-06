"""
Create a python api which accepts a genbank or other tab to construct json for cgview.js

"""
import json
import random
import string
from Bio import SeqIO
from os.path import dirname

template_file = f"{dirname(__file__)}/cgview_template.json"

def get_random():
    return "".join(
        random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase)
        for i in range(4)
    )


default_colors = [
    "rgb(153,153,153)",
    "rgb(247,129,191)",
    "rgb(166,86,40)",
    "rgb(255,255,51)",
    "rgb(255,127,0)",
    "rgb(152,78,163)",
    "rgb(77,175,74)",
    "rgb(55,126,184)",
    "rgb(228,26,28)",
]


class cgview(object):
    def __init__(self,name) -> None:
        self.data = json.load(open(template_file))
        self.legend2color = dict()
        self.set_name(name)
    def set_name(self, name):
        self.data["cgview"]["name"] = name
        self.data["cgview"]["captions"][0]["name"] = name
        
    def set_legend(self, name,color=None):
        if name not in self.legend2color:
            self.legend2color[name] = self.color_palette()
            if color is None:
                color = self.color_palette(name)
            _d = {"name": name, "swatchColor": str(color), "decoration": "arc"}
            self.data["cgview"]["legend"]["items"].append(_d.copy())
        else:
            return

    def color_palette(self, name=None):
        if name is not None:
            return self.legend2color[name]
        else:
            c = list(set(default_colors).difference(set(self.legend2color.values())))
            if len(c) == 0:
                print("no color left")
            return random.choice(c)
        
    def check_fea(self,in_dict):
        in_dict = in_dict.copy()
        if not in_dict.get('contig'):
            raise IOError('No contig information')
        if not in_dict.get('strand'):
            in_dict['strand'] = 1
        if not in_dict.get('favorite'):
            in_dict['favorite'] = False
        return in_dict 
    
    def set_feature(self, info_list, legend_title):
        unique_track_num = legend_title + "_" + get_random()
        keys = [
            "name",
            "type",
            "start",
            "stop",
            "strand",
            "source",
            "legend",
            "contig",
            "tags",'favorite'
        ]
        for info in info_list:
            _d = {}
            for k in keys:
                _d[k] = info.get(k, "")

            _d["legend"] = legend_title
            _d["source"] = unique_track_num
            self.set_legend(legend_title)
            
            _d = self.check_fea(_d)
            self.data["cgview"]["features"].append(_d)
        self.set_tracks(unique_track_num, legend_title)
        
    def set_contig(self,in_fna):
        _d = {"name": "",
              "orientation": "+",
              "length": 0,
              "seq": ""
                    }
        for r in SeqIO.parse(in_fna,'fasta'):
            _d = {"name": r.id,
              "orientation": "+",
              "length": len(r.seq),
              "seq": str(r.seq)
                    }
            self.data["cgview"]["sequence"]['contigs'].append(_d.copy())
            
    def set_tracks(self, unique_track_num, name):
        _d = {
            "name": name,
            "separateFeaturesBy": "none",
            "position": "outside",
            "thicknessRatio": 1,
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": unique_track_num,
        }
        self.data["cgview"]["tracks"].append(_d.copy())

    def set_others(self,):
        pass
    
    def write(self,outfile):
        with open(outfile,'w') as f1:
            json.dump(self.data,f1)

if __name__ == '__main__':
    vb = cgview('B4')
    vb.set_contig('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/B4/09_prokka/B4.fna')
    vb.set_legend('IS','#C62828')
    vb.set_legend('prophage','#3D5AFE')
    vb.set_legend('denitrification','#2E7D32')
    
    import pandas as pd
    _df = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/IS_info.tsv',sep='\t')
    _df.loc[:,'name'] = _df['cluster']
    _df.loc[:,'type'] = 'IS'
    subdf = _df.loc[_df['genome']=='B4',:]
    m = {'+':1,'-':-1}
    subdf.loc[:,'strand'] = [m.get(_,1) for _ in subdf.loc[:,'strand']]
    d = list(subdf.to_dict(orient='index').values())
    vb.set_feature(d,'IS')
    T = 'B4'
    _df = pd.read_csv('/home-user/thliao/project/coral_ruegeria/data_processing/annotations/phage_regions.tab',sep='\t')
    subdf = _df.loc[_df['Genome']==T,:]
    subdf.loc[:,'type'] = 'prophage'
    subdf.loc[:,'stop'] = list(subdf['End'])
    subdf.loc[:,'start'] = list(subdf['Start'])
    subdf.loc[:,'contig'] = [f"{T}_{_}" for _ in subdf['Contig']]
    d = list(subdf.to_dict(orient='index').values())
    vb.set_feature(d,'prophage')
    _df = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/cgviewer/B4_deni.tsv')
    _df.loc[:,'favorite'] = True
    d = list(_df.to_dict(orient='index').values())
    vb.set_feature(d,'denitrification')
    vb.write('./B4.json')
    
    