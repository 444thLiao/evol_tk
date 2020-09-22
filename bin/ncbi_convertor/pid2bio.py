"""
This script is mainly for query requested protein accession id and retrieve relative BioSample and BioProject.
Of course including its relative metadata.

"""
import os
from os.path import exists, dirname

import click
import pandas as pd
from tqdm import tqdm

from bin.ncbi_convertor import parse_id, NCBI_convertor, get_bioproject, get_biosample


def genomeID2Bio(genome_IDs):
    """
    bio include BioSample and corresponding BioProject
    """
    tqdm.write("with genome ID, start to retrieve genome information")
    convertor = NCBI_convertor(genome_IDs, db='assembly')
    convertor.get_GI()
    convertor.get_db_summary()
    aid2info = convertor.dbsummary

    bs_list = list(set([_.get('BioSampleAccn')
                        for _ in aid2info.values()
                        if _.get('BioSampleAccn')]))
    bp_list = list(set([_.get('BioprojectAccn')
                        for _ in aid2info.values()
                        if _.get('BioprojectAccn')]))
    tqdm.write("retrieving relative Bioproject and its Biosample info")
    bp2info = get_bioproject(bp_list)
    bs2info = get_biosample(bs_list)
    return aid2info, bp2info, bs2info


def main(infile, ofile, db='protein', force=False, redo=False):
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    order_id_list, id2annotate = parse_id(infile)

    if db == 'protein':
        convertor = NCBI_convertor(order_id_list,db='protein')
        pid2assembly_dict = convertor.pid2assembly()
        genome_IDs = list([_dict['assembly']
                           for pid,_dict in pid2assembly_dict.items()])
        genome_IDs = [_ for _ in genome_IDs if _]

    elif db == 'genome':
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

        if exists(ofile) and not force:
            tqdm.write("detect existing " + ofile +
                       ' no force param input, so it quit instead of writing.')
            return

        full_df.to_csv(ofile, sep='\t', index=1, index_label='AssemblyAccn raw')
        tqdm.write('finish writing into ' + ofile + ' with tab separator format.')
        
    else:
        raise Exception('too much genomes to process')




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
