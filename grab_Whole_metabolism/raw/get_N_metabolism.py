from bioservices.kegg import KEGG
from collections import defaultdict
k = KEGG()

N_pathway_id = "ko00910"

data = k.get(N_pathway_id)
dict_data = k.parse(data)

# get module
list_module = dict_data['MODULE']
module_dict = defaultdict(dict)
for m in list_module:
    module_data = k.parse(k.get(m))
    module_dict[m]['data'] = module_data
    module_dict[m]['metadata'] = list_module[m]

# get orthology/enzyme
# why don't we use enzyme
# problematic case:
# K00372+K00360
# assimilatory nitrate reductase [EC:1.7.99.-] [RN:R00798]
m2orthology = defaultdict(dict)
for m,mdata in module_dict.items():
    _cache_dict = mdata['data']
    names = [_.split('+') for _ in _cache_dict['ORTHOLOGY']]
    m2orthology[m]['metadata'] = _cache_dict['ORTHOLOGY']
    for each_names in names:
        m2orthology[m]['orthology'] = {}
        for name in each_names:
            m2orthology[m]['orthology'][name] = k.parse(k.get(name))

#
# res = k.parse_kgml_pathway(N_pathway_id)
#
# orthology = [_ for _ in res['entries'] if _['type'] == 'ortholog']