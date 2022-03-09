"""
Implement a python usable library for API provided by NCBI recently
source url : https://api.ncbi.nlm.nih.gov/datasets/v1alpha/

Up to date, it doesn't show any special features need to be used. 

initiated by Tianhua liao
at 20200703
"""
import requests

base_url = "https://api.ncbi.nlm.nih.gov/datasets/v1alpha"
r = requests.get(base_url)


def get_assembly_des_by_tax(tax_name):
    "Assembly descriptions by taxonomic name (scientific or common name at any tax rank)"
    new_url = base_url + '/assembly_descriptors/organism/' + tax_name
    r = requests.get(new_url)
    data = r.json()
    num_return = data['total_count']
    datasets = data['datasets']
    pass


def get_assembly_des_by_ass(tax_name):
    "Assembly descriptions by taxonomic name (scientific or common name at any tax rank)"
    new_url = base_url + '/assembly_descriptors/organism/' + tax_name
    r = requests.get(new_url)
    data = r.json()
    num_return = data['total_count']
    datasets = data['datasets']
    pass
