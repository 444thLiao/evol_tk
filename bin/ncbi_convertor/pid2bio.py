"""
This script is mainly for query requested protein accession id and retrieve relative BioSample and BioProject.
Of course including its relative metadata.

"""
import io
import os
from os.path import exists, dirname

import click
import pandas as pd
from Bio import Entrez
from tqdm import tqdm

from bin.ncbi_convert import edl, parse_id
from bin.ncbi_convert.pid2GI import pid2GI
from bin.ncbi_convert.pid2genome import pid2genome_assembly
from global_search.thirty_party.metadata_parser import parse_bioproject_xml, parse_biosample_xml, parse_assembly_xml


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
    _results = []
    while failed:
        _results, failed = edl.esearch(db='assembly',
                                       ids=genome_IDs,
                                       result_func=lambda x: Entrez.read(
                                           io.StringIO(x))['IdList'],
                                       batch_size=1)
    all_assembly_GI = results[::] + _results
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

    return gid2assembly_info, bp2info, bs2info


def main(infile, ofile, start_at='protein', force=False, redo=False):
    if not exists(dirname(ofile)) and dirname(ofile):
        os.makedirs(dirname(ofile))

    order_id_list, id2annotate = parse_id(infile)
    if start_at == 'protein':
        id2gi = {}
        if isinstance(id2annotate[order_id_list[0]], dict):
            # it is a dict, so it contains other infomation or implemented GI. it may be passed over.
            if 'GI' in id2annotate[order_id_list[0]]:
                print("provided file already contains `GI` column(doesn't check the validation/completeness). Giving `force` param to overwrite/implement it. ")
                if not force:
                    id2gi = {k: id2annotate[k]['GI'] for k in order_id_list}
            # todo: re-implemented original infomation into `ofile` from `infile`
        else:
            # no header, just a list of IDs
            pass
        if not id2gi:
            id2gi = pid2GI(order_id_list, redo=redo)
        pid2assembly_dict = pid2genome_assembly(id2gi, redo=redo)
        with_genome_pid2ass = {k: v for k, v in pid2assembly_dict.items() if set(v) != {
            ''}}
        genome_IDs = set([_ for v in with_genome_pid2ass.values() for _ in v])
        genome_IDs = [_ for _ in genome_IDs if _]
    elif start_at == 'genome':
        genome_IDs = order_id_list[::]
    else:
        raise SyntaxError('wrong input of start_at')
    gid2assembly_info, bp2info, bs2info = genomeID2Bio(genome_IDs)

    # if not too big, use panda to concat them
    # else...... too complicated...pass it
    if len(gid2assembly_info) <= 15000:
        ginfo_df = pd.DataFrame.from_dict(gid2assembly_info, orient='index')
        # ginfo_df.index = ginfo_df.iloc[:,0]
        bp_df = pd.DataFrame.from_dict(bp2info, orient='index')
        bs_df = pd.DataFrame.from_dict(bs2info, orient='index')
        _df1 = bp_df.reindex(ginfo_df.loc[:, 'BioprojectAccn'])
        _df1.index = ginfo_df.index
        _df2 = bs_df.reindex(ginfo_df.loc[:, 'BioSampleAccn'])
        _df2.index = ginfo_df.index
        full_df = pd.concat([ginfo_df,
                             _df1,
                             _df2], axis=1)
        full_df = full_df.applymap(lambda x: x.replace('\n', ' ')
        if isinstance(x, str) else x)
        full_df = full_df.drop(['GI', 'relative biosample'], axis=1)
    else:
        raise Exception('too complicated to deal with it')

    if exists(ofile) and not force:
        tqdm.write("detect existing " + ofile +
                   ' no force param input, so it quit instead of writing.')
        return

    full_df.to_csv(ofile, sep='\t', index=1, index_label='AssemblyAccn raw')
    tqdm.write('finish writing into ' + ofile + ' with tab separator format.')


@click.command()
@click.option('-i', 'infile', help='input file which contains protein accession id ')
@click.option('-o', 'ofile', help='output file')
@click.option('-s', 'start_at', help='start from `protein` or `genome` ID.  etc, protein id maybe like `CBH97221.1`. genome ID should like `GCF_900176205.1` ', default='protein',
              required=False)
@click.option('-f', 'force', help='force overwrite?', default=False, required=False, is_flag=True)
@click.option('-redo', 'redo', help='use cache or not? default is use the cache.', default=False, required=False, is_flag=True)
def cli(infile, ofile, force, redo, start_at):
    start_at = start_at.strip().lower()
    if start_at not in ['genome', 'protein']:
        raise IOError('Unexpected params of start_at. giving `%s`?? ' % start_at)
    main(infile, ofile, start_at, force, redo)


if __name__ == "__main__":
    cli()
