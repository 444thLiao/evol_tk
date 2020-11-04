from bs4 import BeautifulSoup
from collections import defaultdict


def parse_bioproject_xml(xml_text):
    result_bucket = []
    soup = BeautifulSoup(xml_text, 'xml')
    split_out = soup.find_all("DocumentSummary")
    for each_record in split_out:
        bioproject2info = defaultdict(dict)
        body = each_record.Project
        uid = each_record['uid']
        if body is None:
            result_bucket.append(uid)
            continue
        major_parts = {_.name: _ for _ in body.children if _.name}
        PID_part = major_parts['ProjectID']
        bioproject_id = PID_part.ArchiveID['accession']
        bioproject_gid = int(PID_part.ArchiveID['id'])
        PID_des = major_parts['ProjectDescr']
        des_text = PID_des.Description.text if PID_des.Description is not None else ''
        des_title = PID_des.Title.text if PID_des.Title is not None else ''
        PID_pubmed_id = ';'.join([_['id']
                                  for _ in PID_des.find_all('Publication')])
        list_biosample = PID_des.find_all('LocusTagPrefix')
        biosample_text = ';'.join([_.get('biosample_id', 'no biosample ID')
                                   for _ in list_biosample])

        PID_type = major_parts['ProjectType']
        biological_properties = PID_type.find('BiologicalProperties')
        if biological_properties is not None:
            env_data = biological_properties.find('Environment')
            if env_data:
                for _data in env_data:
                    if _data != '\n':
                        bioproject2info[bioproject_id][_data.name] = _data.text
            Morphology_data = biological_properties.find('Morphology')
            if Morphology_data:
                for _data in Morphology_data:
                    if _data != '\n':
                        bioproject2info[bioproject_id][_data.name] = _data.text
        target_pro = PID_type.find('target')
        if target_pro:
            for _data,v in target_pro.attrs.items():
                bioproject2info[bioproject_id][_data] = v
        bioproject2info[bioproject_id]['GI'] = bioproject_gid
        bioproject2info[bioproject_id]['description'] = des_text
        bioproject2info[bioproject_id]['title'] = des_title
        bioproject2info[bioproject_id]['pubmed GI'] = PID_pubmed_id
        bioproject2info[bioproject_id]['relative biosample'] = biosample_text
        bioproject2info[bioproject_id]['number of biosamples'] = len(
            list_biosample)
        pl = each_record.ProjectLinks
        if pl is not None:
            memids_DOM = [_ for _ in pl.find_all("MemberID")]
            memberids = ';'.join([_['id'] for _ in memids_DOM])
            member_accessions = ';'.join([_['accession'] for _ in memids_DOM])
            bioproject2info[bioproject_id]['member bioproject GI'] = memberids
            bioproject2info[bioproject_id]['member bioproject accession'] = member_accessions
        result_bucket.append(bioproject2info)
    return result_bucket


def parse_biosample_xml(xml_text):
    result_bucket = []
    soup = BeautifulSoup(xml_text, 'xml')
    split_out = soup.find_all("BioSample")
    for each_record in split_out:
        biosample2info = defaultdict(dict)
        uid = each_record['id']
        accession = each_record['accession']
        biosample2info[accession]['access status'] = each_record.get('access','')
        biosample2info[accession]['last_update'] = each_record.get('last_update','')
        biosample2info[accession]['publication_date'] = each_record.get('publication_date','')
        biosample2info[accession]['submission_date'] = each_record.get('submission_date','')

        for dom in each_record.find_all('Id'):
            for attr in dom.attrs:
                if attr.startswith('db'):
                    biosample2info[accession][dom[attr]] = dom.text
        for dom in each_record.find('Description').children:
            if dom.name is not None:
                biosample2info[accession][dom.name] = dom.text
        for dom in each_record.find_all('Attribute'):
            biosample2info[accession]['attribute:' +
                                      dom['attribute_name']] = dom.text
        result_bucket.append(biosample2info)
    return result_bucket


def parse_assembly_xml(xml_text):
    result_bucket = []
    soup = BeautifulSoup(xml_text, 'xml')
    split_out = soup.find_all("DocumentSummary")
    for each_record in split_out:
        assembly2info = defaultdict(dict)
        uid = each_record['id']
        aid = each_record.find('AssemblyAccession').text.strip().strip('\n')
        info_get_ = ['Genbank',
                     'SpeciesName',
                     'Isolate',
                     "ExclFromRefSeq",
                     "Coverage",
                     'Infraspecie',
                     'SpeciesTaxid',
                     'AssemblyStatus',
                     'BioSampleAccn',
                     'BioSampleId',
                     'FtpPath_GenBank',
                     'FtpPath_RefSeq']
        assembly2info[aid]['GI'] = uid
        for key in info_get_:
            _cache = each_record.find(key)
            if _cache:
                if key == 'Infraspecie':
                    _cache = _cache.find('Sub_value')
                    key = 'Isolate'
                    if _cache is None:
                        continue
                    assembly2info[aid][key] = _cache.text.strip().strip('\n')
                elif key == 'Genbank':
                    key = 'AssemblyAccession'
                    assembly2info[aid][key] = _cache.text.strip().strip('\n')
                else:
                    assembly2info[aid][key] = _cache.text.strip().strip('\n')
        
        bioproject_total = each_record.find_all('GB_BioProjects')
        if bioproject_total:
            bioproject_total = bioproject_total[0]
            if bioproject_total.find('BioprojectAccn'):
                bp_accn = bioproject_total.find('BioprojectAccn').text.strip().strip('\n')
                assembly2info[aid]['BioprojectAccn'] = bp_accn
            if bioproject_total.find('BioprojectId'):
                bp_gi = bioproject_total.find('BioprojectId').text.strip().strip('\n')
                assembly2info[aid]['BioprojectId'] = bp_gi
        result_bucket.append(assembly2info)
    return result_bucket
        
        