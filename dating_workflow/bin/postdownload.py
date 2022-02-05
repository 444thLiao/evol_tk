"""
It mainly for running prokka after download genomes.
It includes
1. check the missing ID
2. perform prokka and rename the ID
3. generate metadata.csv
"""
import multiprocessing as mp
import os
from glob import glob
from os.path import join, exists, basename, dirname
from subprocess import check_call

import click
from tqdm import tqdm

from dating_workflow.step_script import convert_genome_ID


def run_cmd(cmd):
    check_call(cmd, shell=True,
               stderr=open('/dev/null', 'w'),
               stdout=open('/dev/null', 'w'))


def get_faa_from_prokka_r(infile,
                          odir,
                          sample_name,
                          prokka_p,
                          return_cmd=False,
                          thread_per_prokka=0):
    locustag = convert_genome_ID(sample_name)
    oprefix = f"{odir}/{sample_name}/{sample_name}"
    if exists(f'{oprefix}.faa'):
        # To existing prokka output, it would not rerun it.
        return f'{oprefix}.faa'

    cmd = f'{prokka_p} {infile} --outdir {odir}/{sample_name} --force --prefix {sample_name} --locustag {locustag} --cpus {thread_per_prokka} '
    if return_cmd:
        return cmd
    run_cmd(cmd)
    return f'{oprefix}.faa'


def cli(indir,
        odir=None, 
        tmp_dir=None,
        reformatted_name=True,
        force=False,
        prokka_p=" `which prokka`",
        thread=5,
        all_ids=None,
thread_per_prokka=0,
dry_run=False,
protein=False
        ):
    """
    It would use the downloaded protein first.
    If it doesn't exist, it will perform prokka to predict genes.
    :param indir: ./genbank
    :param odir: ./modified_data/genome_protein_files/
    :param tmp_dir: ./modified_data/prokka_o/
    :param reformatted_name: reformatted the output file or not
    :param force: overlap the output file or not
    :param force_prokka: deprecated
    :param prokka_p: the exec path of prokka. Default using `which prokka` to retrieve it path
    :return:
    """

    # process the directory
    if odir is None:
        odir = './genome_protein_files'
    if tmp_dir is None:
        tmp_dir = join(odir, 'tmp')
        # tmp_dir = join(odir, 'tmp')
    if not exists(tmp_dir):
        os.makedirs(tmp_dir, exist_ok=True)

    tqdm.write(f'iterating the {indir}')
    
    # all_dir = [_
    #     for _ in tqdm(glob(join(indir, '**', 'GC*', '*.fna.gz'))
    #                     )
    #     ]
    if all_ids:
        paths = os.listdir(indir)
        all_dir = []
        for aid in all_ids:
            for p in paths:
                p = indir+'/'+p
                if exists(join(p,aid)):
                    d = glob(join(p, aid, '*.fna.gz'))
                    if d:
                        all_dir.append(d[0])
                    break
    else:
        all_dir = [_
            for _ in tqdm(glob(join(indir, '**', 'GC*', '*.fna.gz'))
                            )
            ]        
        all_dir = [_ for _ in all_dir if basename(dirname(_)) in all_ids]
        found = [basename(dirname(_)) for _ in all_dir]
        not_found = list(set(all_ids).difference(set(found)))
        if not_found:
            tqdm.write(f"{len(not_found)} are not found!")
            
    tqdm.write("gunzip fna file and collect jobs")
    jobs = []
    jobs2 = []
    for p_dir in tqdm(all_dir):
        p_dir = dirname(p_dir)
        p_files = glob(join(p_dir, '*.faa.gz'))
        if protein and p_files:
            # if protein is given, it indicate that it use the downloaded protein file directly instead of the protein file annotated by prokka
            faa = p_files[0]
            run_cmd(f"gunzip -d -c {faa} > {faa.replace('.gz', '')}")
            ofile = join(odir, basename(p_dir)) + '.faa'
            cmd = f"ln -s `realpath {faa.replace('.gz', '')}` {ofile}"
            jobs.append(cmd)
            continue
        
        ofile = join(odir, basename(p_dir)) + '.faa'
        if exists(ofile) and not force:
            # if the output faa exists and not force, pass  it
            continue
        # if not p_files:
        # it haven't protein files
        # use prokka to predice gene
        fna_file = glob(join(p_dir, '*.fna.gz'))[0]
        new_fna = fna_file.replace('.gz', '')
        if not exists(new_fna):
            run_cmd(f'gunzip -d -c {fna_file} > {new_fna}')
        sample_name = basename(dirname(fna_file))
        prokka_cmd = get_faa_from_prokka_r(infile=new_fna,
                                           odir=tmp_dir,
                                           sample_name=sample_name,
                                           prokka_p=prokka_p,
                                           return_cmd=True,
                                           thread_per_prokka=thread_per_prokka,
                                           )
        if exists(prokka_cmd):
            # output is a file instead of cmd which means it has prokka_output already.
            prokka_ofile = prokka_cmd
            jobs2.append(f"ln -s `realpath {prokka_ofile}` {ofile}")
            continue
        else:
            jobs.append(prokka_cmd)

        # collect prokka output file
        prokka_ofile = f"{tmp_dir}/{sample_name}/{sample_name}.faa"
        jobs2.append(f'ln -s `realpath {prokka_ofile}` {ofile}')

        # if p_file.endswith('.gz') and exists(prokka_ofile):
        #     run_cmd(f'gunzip -d -c {p_file} >{ofile}')
        # elif  exists(prokka_ofile):
        #     run_cmd(f'cat {p_file} >{ofile}')
        # else:
        #     # print(p_file,ofile)
        #     pass
    if dry_run:
        with open('./cmds','w') as f1:
            f1.write('\n'.join(jobs))
        return
    else:
        tqdm.write('run prokka')
        with mp.Pool(processes=thread) as tp:
            r = list(tqdm(tp.imap(run_cmd, jobs), total=len(jobs)))

    tqdm.write('run soft link to save space')
    for j in tqdm(jobs2):
        run_cmd(j)
    # processed_ids = [_.split('/')[-1].strip().replace('.faa', '') for _ in jobs2]
    # format protein id
    # one thing need to be noted. The ID produced by prokka would not directly equal to the formatted name here. Because the number of prokka would include the rRNA gene. the formatted name here would not.
    
    # dangerous behaviour
    # if reformatted_name:
    #     name_map = {}
    #     for p in tqdm(glob(join(odir, '*.faa'))):
    #         name = basename(p).replace('.faa', '')
    #         if not force and name not in processed_ids:
    #             continue
    #         locus_prefix = convert_genome_ID(name)
    #         records = []
    #         for idx, record in enumerate(SeqIO.parse(p, format='fasta')):
    #             new_name = locus_prefix + '_{:0>5}'.format(idx + 1)
    #             if record.id.split('_')[0] == locus_prefix:
    #                 continue
    #             name_map[record.id] = new_name
    #             record.id = new_name
    #             records.append(record)
    #         SeqIO.write(records,
    #                     open(p, 'w'),
    #                     format='fasta-2line')


@click.command(help="""
It mainly for preprocess downloaded genomes using prokka or its original proteins files.
It includes
1. check the missing ID
2. perform prokka and rename the ID
3. generate metadata.csv
""")
@click.option("-i", "indir", help="input directory [./genbank]", default="./genbank", )
@click.option("-o", "odir", help="input directory [./genome_protein_files]", default="./genome_protein_files")
@click.option("-tmp", "tmp_dir", help='For saving time and space, you could assign tmp_dir [./tmp]',
              default=None)
@click.option("-gl", "genome_list", default=None,
              help="It will read 'selected_genomes.txt', please prepare the file, or indicate the alternative name or path. It could be None. If you provided, you could use it to subset the aln sequences by indicate names.")
@click.option('-f', 'force', help='overwrite? mainly for prokka', default=False, required=False, is_flag=True)
@click.option('-protein', 'protein', help='unzip protein file and link it to odir', 
default=False, required=False, is_flag=True)
@click.option('-np', 'num_parellel', default=5, help="num of processes could be parellel.. default is 10")
@click.option('-t', 'thread_per_prokka', default=0, help="num of threads prokka used. default 0 means all threads of a server.")
@click.option('-dry_run', 'dry_run', default=False, is_flag=True, help="output commands to './cmds' needed to run. It only includes the prokka part. You need to soft link/move them to the odir don't run it. ")
def main(indir, odir, tmp_dir, genome_list, num_parellel, force,thread_per_prokka,dry_run,protein):
    if not exists(indir):
        raise IOError("input dir doesn't exist")
    if not exists(odir):
        os.makedirs(odir, exist_ok=True)

    ## The metadata would be parsed at the downloading step.
    # for indir in ['./genbank', './refseq']:
    #     if exists(indir):
    # name = indir.split('/')[-1]
    # base_tab = expanduser(f'~/.cache/ncbi-genome-download/{name}_bacteria_assembly_summary.txt')
    # all_g_ids = set([basename(_)
    #                  for _ in glob(join(indir, 'bacteria', '*'))])
    # # from downloaded dir
    if genome_list is None:
        all_ids = None
    else:
        all_ids = open(genome_list).read().split('\n')
        all_ids = [_ for _ in all_ids if _]
        all_ids = set(all_ids)
    # # from id list
    # metadatas = open(base_tab).read().split('\n')
    # rows = [_
    #         for _ in tqdm(metadatas)
    #         if _.split('\t')[0] in all_g_ids]
    #
    # f1 = open(join(odir,'metadata.csv'),'w')
    # f1.write(metadatas[1].strip('# ') + '\n')
    # f1.write('\n'.join(rows))
    # f1.close()
    # if set(all_ids) != set(all_g_ids):
    #     print('inconsistent id, missing ids: ' + '\n'.join(all_ids.difference(all_g_ids)))

    cli(indir, odir,
        tmp_dir=tmp_dir,
        force=force,
        thread=num_parellel,
        all_ids=all_ids,
        thread_per_prokka=thread_per_prokka,
        protein=protein,
        dry_run=dry_run)


if __name__ == "__main__":
    main()
    # python3 /home-user/thliao/script/evol_tk/dating_workflow/toolkit/postdownload.py -i ./genbank -o ./genome_protein_files
