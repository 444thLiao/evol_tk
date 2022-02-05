"""
Try to implement a API for 'https://api.ncbi.nlm.nih.gov/datasets/v1alpha/'
startted at 20220205
"""


r['assembly']['assembly_accession']
r['assembly']['paired_assembly_accession']
gids = [r['assembly'].get('paired_assembly_accession','') for r in b['assemblies']]
gids2 = [r['assembly'].get('assembly_accession','') for r in b['assemblies']]
gids = [_['']]