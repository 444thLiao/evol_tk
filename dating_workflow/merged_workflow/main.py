from os.path import *
from subprocess import check_call

script_dir = "/home-user/thliao/script/evolution_relative/dating_workflow"


def run(cmd):
    check_call(cmd,
               shell=True)


def prepare_dating(genome_data_dir, genome_list, ofile):
    gids = open(genome_list).read().split('\n')
    gids = [_ for _ in gids if _]
    num_gids = len(set(gids))

    seq_odir = "./cog25_single/seq"
    cog25_first = f"python3 {script_dir}/step_script/extract_cog25.py '{genome_data_dir}/*.faa' ./cog25_annotate {seq_odir}"

    alignment_odir = f"./cog25_single/{num_gids}g_aln"
    alignment_cog25 = f"python3 {script_dir}/bin/batch_run/batch_mafft.py -i {seq_odir} -s faa -o {alignment_odir} -f -m ginsi -gl {genome_list}"
    trimal_cog25 = f"python3 {script_dir}/bin/batch_run/batch_trimal.py -i {alignment_odir} -o {alignment_odir}"
    concat_cog25 = f"python3 {script_dir}/toolkit/concat_aln.py -i {alignment_odir} -ct phy -gl {genome_list} -o {ofile} -s trimal -no_graph"

    if not exists(seq_odir):
        run(cog25_first)

    run(alignment_cog25)
    run(trimal_cog25)
    run(concat_cog25)


if __name__ == "__main__":
    prepare_dating("./rawdata/genome_protein_files",
                   "./trees/iqtree/dating_g_from_over80p_final.list",
                   "./dating_for/phy_files/168g_concat.trimal")
