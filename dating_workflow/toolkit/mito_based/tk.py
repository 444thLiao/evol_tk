from collections import defaultdict
from glob import glob
import os
from os.path import *
from Bio import SeqIO
from tqdm import tqdm
import pandas as pd
from dating_workflow.step_script.quick_sampling import *
from bin.batch_run.batch_mafft import remove_3rd
import copy



def get_faa(aid):
    in_p_dir = [
        "/mnt/maple/thliao/data/NCBI/modified_data/direct_protein_files",
        "/mnt/home-backup/thliao/NCBI/euk_db/proteins_sswang",
    ]
    in_p = [glob(f"{_}/{aid}*")[0] for _ in in_p_dir if len(glob(f"{_}/{aid}*")) == 1]
    if in_p:
        in_p = in_p[0]
    else:
        raise IOError(f'{aid} was not found')
    return in_p

def get_ffn(aid):
    faa = f"/mnt/maple/thliao/data/NCBI/modified_data/direct_protein_files/{aid}.faa"
    if exists(faa):
        ffn = realpath(faa).replace('.faa','.ffn')
    else:
        ffn = f"/mnt/home-backup/thliao/NCBI/euk_db/genes_sswang/{aid}.gene"
    if not exists(ffn):
        raise IOError(f'{aid} was not found')
    return ffn


def anno_repr(all_ids, _outdir):
    db = "/mnt/maple/thliao/data/protein_db/mito-based_orthologs/protein/representative.protein"
    cmds = []
    # _outdir = f'./dating/mito/anno/{dataset}/anno_tab'
    if not exists(_outdir):
        os.makedirs(_outdir)
    for aid in tqdm(all_ids):
        in_p = get_faa(aid)
        ofile = f"{_outdir}/{aid}__{basename(db)}.tab"
        cmd = f"blastp -query {in_p} -db {db} -outfmt 6 -out {ofile} -evalue 1e-10 -max_hsps 1 -num_threads 20"
        if not exists(ofile) or getsize(ofile) == 0:
            cmds.append(cmd)
    return cmds


def parse_anno(_outdir, locus2len, name):
    """filtering critera 
    
    evalue < 1e-10 and Query coverage >= 80%
    
    """
    mi2genome2list_locus = defaultdict(dict)
    mi2genome2locus = defaultdict(dict)
    genome2mito = defaultdict(list)
    for tab in tqdm(glob(join(_outdir, "*__*.tab"))):
        aid = basename(tab).split("__")[0]
        if os.path.getsize(tab) != 0:
            _df = pd.read_csv(tab, sep="\t", header=None)
            _df = _df.groupby(0).head(1)  # make each locus only be mapped to a gene
            _df.loc[:, "gene"] = [_ for _ in _df[1]]
            _df.loc[:, "align len"] = abs(_df[7] - _df[6])
            _df.loc[:, "total len"] = [locus2len.get(_, 0) for _ in _df[0]]
            _df.loc[:, "cov"] = _df.loc[:, "align len"] / _df.loc[:, "total len"] * 100
            _df = _df.loc[
                (_df[10] <= 1e-10) & (_df["cov"] >= 80), :
            ]  # add coverage over 80
            _df = _df.sort_values(10)  # ascending
            gb = _df.groupby("gene")
            for gene, locus_list in gb.groups.items():
                locus_list = list(_df.loc[locus_list, 0])
                all_locus = list(set(locus_list))
                mi2genome2list_locus[gene][aid] = ";".join(all_locus)
            sub_df = gb.head(1)
            gene2locus = dict(zip(sub_df["gene"], sub_df[0]))
            for gene, locus in gene2locus.items():
                mi2genome2locus[gene][aid] = locus
                genome2mito[aid].append(gene)
    _df = pd.DataFrame.from_dict(mi2genome2list_locus).fillna("-")
    _bin_df = _df.applymap(lambda x: len([_ for _ in x.split(";") if _ != "-"]))
    with pd.ExcelWriter(join(_outdir, "..", f"{name}.xlsx")) as f1:
        _df.to_excel(f1, "locus info", index=1, index_label="Assembly IDs")
        _bin_df.to_excel(f1, "number of locus", index=1, index_label="Assembly IDs")
    return genome2mito, mi2genome2locus


def output_seqs(odir, genome2gene, gene2genome2locus, mode="prot"):
    if not exists(odir):
        os.makedirs(odir)
    all_genomes = set()
    gene2seqs = defaultdict(list)
    locus2seq = {}
    for aid in genome2gene:
        if mode == "prot":
            locus2seq.update(
                {r.id: r for r in list(SeqIO.parse(get_faa(aid), "fasta"))}
            )
        elif mode == "nucl":
            locus2seq.update(
                {r.id: r for r in list(SeqIO.parse(get_ffn(aid), "fasta"))}
            )
    for gene, _d in tqdm(gene2genome2locus.items()):
        for genome, locus in _d.items():
            r = locus2seq[locus]
            r = copy.deepcopy(r)
            r.id = r.name = r.description = ""
            r.id = genome
            gene2seqs[gene].append(r)
            all_genomes.add(genome)
    with open(odir + "/used_genomes.list", "w") as f1:
        f1.write("\n".join(list(all_genomes)))
    for mi, seqs in gene2seqs.items():
        with open(join(odir, f"{mi}.faa"), "w") as f1:
            if mode == "nucl":
                SeqIO.write([remove_3rd(_) for _ in seqs], f1, "fasta-2line")
            elif mode == "prot":
                SeqIO.write(seqs, f1, "fasta-2line")
            else:
                tqdm.write("noting happen")


def run_concat(odir,ofile,suffix='.faa'):
    cmds = []
    for faa in glob(join(odir,'*'+suffix)):
        if not exists(faa.replace(suffix,'.aln')):
            cmd = f"ginsi {faa} > {faa.replace(suffix,'.aln')}; trimal -in {faa.replace(suffix,'.aln')} -out {faa.replace(suffix,'.trimal')} -automated1 -keepheader"
            cmds.append(cmd)
    cmds.append(f"""python3 /home-user/thliao/script/evol_tk/dating_workflow/bin/concat_aln.py -i {odir} -o {ofile} -s trimal -simple -gl {odir+'/used_genomes.list'} -ct both -no_graph 
    """)
    return cmds