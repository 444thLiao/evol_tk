"""
This script is mainly for query requested protein accession id and retrieve relative BioSample and BioProject.
Of course including its relative metadata.

"""
from bin.ncbi_convert import edl, access_intermedia, parse_id
from bin.ncbi_convert.pid2GI import pid2GI
from bin.ncbi_convert.pid2tax import GI2tax
from bin.ncbi_convert.pid2genome import pid2genome_assembly
from global_search.thirty_party.metadata_parser import parse_bioproject_xml,parse_biosample_xml,parse_assembly_xml
from os.path import exists, join, dirname
from tqdm import tqdm
from Bio import Entrez
import io
import os
import click
import pandas as pd


def get_bioproject(bp_list):
    results, failed = edl.esearch(db='bioproject',
                                    ids=bp_list,
                                    result_func=lambda x: Entrez.read(io.StringIO(x))['IdList'])
    all_GI = results[::]
    results, failed = edl.efetch(db='bioproject',
                                    ids=all_GI,
                                    retmode='xml',
                                    retype='xml',
                                    result_func=lambda x: parse_bioproject_xml(x))
    bp2info = {}
    for _ in results:
        if isinstance(_, dict):
            bp2info.update(_)
    return bp2info

def get_biosample(bs_list):
    results, failed = edl.esearch(db='biosample',
                                    ids=bs_list,
                                    result_func=lambda x: Entrez.read(io.StringIO(x))['IdList'])
    all_GI = results[::]
    results, failed = edl.efetch(db='biosample',
                                    ids=all_GI,
                                    retmode='xml',
                                    retype='xml',
                                    result_func=lambda x: parse_biosample_xml(x))
    bs2info = {}
    for _ in results:
        if isinstance(_, dict):
            bs2info.update(_)
    return bs2info     
            
def genomeID2Bio(genome_IDs):
    """
    bio include BioSample and corresponding BioProject
    """
    tqdm.write("with genome ID, start to retrieve genome information")
    results, failed = edl.esearch(db='assembly',
                                  ids=genome_IDs,
                                  result_func=lambda x: Entrez.read(io.StringIO(x))['IdList'])
    all_assembly_GI = results[::]
    results, failed = edl.esummary(db='assembly',
                                   ids=all_assembly_GI,
                                   result_func=lambda x: parse_assembly_xml(x))
    gid2assembly_info = {}
    for _ in results:
        gid2assembly_info.update(_)
    bs_list = list(set([_.get('BioSampleAccn') 
               for _ in gid2assembly_info.values() 
               if _.get('BioSampleAccn')]))
    bp_list = list(set([_.get('BioprojectAccn') 
               for _ in gid2assembly_info.values() 
               if _.get('BioprojectAccn')]))
    tqdm.write("retrieving relative Bioproject and its Biosample info")
    bp2info = get_bioproject(bp_list)
    bs2info = get_biosample(bs_list)
    
    return gid2assembly_info,bp2info,bs2info
    

def main(infile, ofile, force=False,redo=False):
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
        
    order_id_list, id2annotate = parse_id(infile)
    id2gi = {}
    if isinstance(id2annotate[order_id_list[0]], dict):
        # it is a dict, so it contains other infomation or implemented GI. it may be passed over.
        if 'GI' in id2annotate[order_id_list[0]]:
            print("provided file already contains `GI` column(doesn't check the validation/completeness). Giving `force` param to overwrite/implement it. ")
            if not force:
                id2gi = {k:id2annotate[k]['GI'] for k in order_id_list}
        # todo: re-implemented original infomation into `ofile` from `infile`
    else:
        # no header, just a list of IDs
        pass
    if not id2gi:
        id2gi = pid2GI(order_id_list,redo=redo)
    pid2assembly_dict = pid2genome_assembly(id2gi,redo=redo)
    with_genome_pid2ass = {k:v for k,v in pid2assembly_dict.items() if set(v)!={''}}
    genome_IDs = set([_ for v in with_genome_pid2ass.values() for _ in v])
    genome_IDs = [_ for _ in genome_IDs if _]
    
    gid2assembly_info,bp2info,bs2info = genomeID2Bio(genome_IDs)
    
    #if not too big
    if len(gid2assembly_info) <= 5000:
        ginfo_df = pd.DataFrame.from_dict(gid2assembly_info,orient='index')
        bp_df = pd.DataFrame.from_dict(bp2info,orient='index')
        bs_df = pd.DataFrame.from_dict(bs2info,orient='index')
        full_df = pd.concat([ginfo_df,
                             bp_df.reindex(ginfo_df.loc[:,'BioProject']),
                             bp_df.reindex(ginfo_df.loc[:,'BioSample'])])
        
    if exists(ofile) and not force:
        tqdm.write("detect existing "+ofile+' no force param input, so it quit instead of writing.')
        return 
    
    with open(ofile, 'w') as f1:
        print('#accession ID\tGI\tassembly_ID\tnuccore ID\tstart\tend\tstrand', file=f1)
        for pid, nuc_dict in pid2assembly_dict.items():
            GI = id2gi[pid]
            for assembly_id,info in nuc_dict.items():
                print(f'{pid}\t{GI}\t{assembly_id}\t' + '\t'.join(map(str,info)), file=f1)
    tqdm.write('finish writing into ' + ofile)


@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id ')
@click.option('-o', 'ofile', help='output file')
@click.option('-f', 'force', help='force overwrite?', default=False, required=False, is_flag=True)
@click.option('-redo', 'redo', help='use cache or not? default is use the cache.', default=False, required=False, is_flag=True)
def cli(infile, ofile, force,redo):
    main(infile, ofile, force,redo)


if __name__ == "__main__":
    cli()
