import xml.etree.ElementTree as ET
from collections import defaultdict


def getfromname(node, name):
    cache = list(node.iterfind(name))
    if not cache:
        return ''
    else:
        return cache[0].text


def xml2dict(xml_data):
    tree = ET.fromstring(xml_data)
    acc2info = defaultdict(dict)
    for child in tree:
        accession = child.attrib['accession']
        acc2info[accession].update(child.attrib)
        for info in child:
            if info.tag in ['TITLE', 'DESCRIPTION', 'NAME', 'ASSEMBLY_LEVEL', 'GENOME_REPRESENTATION']:
                acc2info[accession][info.tag] = info.text
            elif info.tag in ['TAXON']:
                acc2info[accession]['taxon'] = getfromname(info, 'TAXON_ID')
                acc2info[accession]['sci_name'] = getfromname(info, 'SCIENTIFIC_NAME')
            elif info.tag in ['SAMPLE_REF', 'STUDY_REF']:
                _name = info.tag.split('_')[0]
                info = list(info)[0]
                acc2info[accession][f'{_name}_id'] = getfromname(info, 'PRIMARY_ID')
                acc2info[accession][f'{_name}_v2'] = getfromname(info, 'SECONDARY_ID')
            elif info.tag in ['ASSEMBLY_LINKS']:
                """
                <ASSEMBLY_LINKS>
                <ASSEMBLY_LINK>
                  <URL_LINK>
                    <LABEL>WGS_SET_FLATFILE</LABEL>
                    <URL>ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/ca/CABIWI01.dat.gz</URL>
                  </URL_LINK>
                </ASSEMBLY_LINK>
                """
                for link in info:
                    link = list(link)[0]
                    label = getfromname(link, 'LABEL')
                    url = getfromname(link, 'URL')
                    acc2info[accession][f'url of {label}'] = url
    return acc2info


def cli():
    pass


if __name__ == '__main__':
    cli()