from collections import defaultdict
from glob import glob
from os.path import join, exists, dirname
import os

infaa_dir = "/mnt/ivy/thliao/project/coral_ruegeria/data_processing/all_faa"
suffix = ".faa"
db_path = "/mnt/ivy/thliao/project/coral_ruegeria/extra_genes_db/tda_syn/ref.faa"
odir = "/mnt/ivy/thliao/project/coral_ruegeria/extra_genes_db/tda_syn/annotations"

if not exists(odir):
    os.system(f"mkdir -p {odir}")

cmds = []
for faa in glob(join(infaa_dir, "*" + suffix)):
    gid = faa.split("/")[-1].replace(suffix, "")
    cmds.append(
        f"""blastp -query {faa} -db {db_path} -out {odir}/{gid}.blast -max_target_seqs 10000000 -evalue 1e-20 -num_threads 9  -outfmt '6 qaccver saccver staxid pident length mismatch gapopen qstart qend sstart send evalue bitscore' """
    )

from bin.multiple_sbatch import sbatch_all
sbatch_all(cmds, thread_per_tasks=9, batch_size=30, prefix_name="blast")


# from api_tools.IO_for import _parse_blastp
gene2locus = defaultdict(list)
genome2genes = {}
for out in glob(join(odir, f"*.blast")):
    genome = out.split("/")[-1].replace(".blast", "")
    gid2locus = _parse_blastp(out)
    for gene, locus in gid2locus.items():
        gene2locus[gene].extend(locus)
    genome2genes[genome] = list(gid2locus)


from api_tools.itol_func import to_binary_shape

text = to_binary_shape(genome2genes, same_color="", unfilled_other=True)
with open(dirname(odir) + "/tda_binary.txt", "w") as f1:
    f1.write(text)

genome2genes = defaultdict(list)
map_d = dict(
    zip(
        ["K17486", "K20034", "K20035", "K20036"],
        [
            "dmdA",
            "dmdB",
            "dmdC",
            "dmdD",
        ],
    )
)
sub_d = kegg_df.loc[:, ["K17486", "K20034", "K20035", "K20036"]].to_dict()
for gene, _d in sub_d.items():
    for genome, v in _d.items():
        gene_name = map_d[gene]
        if not pd.isna(v) or not v:
            genome2genes[genome].append(gene_name)
from api_tools.itol_func import to_binary_shape

text = to_binary_shape(genome2genes, same_color="", unfilled_other=True)
with open(dirname(odir) + "/dmd_binary.txt", "w") as f1:
    f1.write(text)
genome2genes = {
    k: ["dmdD"]
    for k, v in kegg_df.loc[:, "K20036"].to_dict().items()
    if not pd.isna(v) or not v
}
from api_tools.itol_func import to_binary_shape

text = to_binary_shape(genome2genes, same_color="", unfilled_other=True)
with open(dirname(odir) + "/dmdD_binary.txt", "w") as f1:
    f1.write(text)
genome2genes = {
    k: ["dddL"]
    for k, v in kegg_df.loc[:, "K16953"].to_dict().items()
    if not pd.isna(v) or not v
}
from api_tools.itol_func import to_binary_shape

text = to_binary_shape(genome2genes, same_color="", unfilled_other=True)
with open(dirname(odir) + "/dddL_binary.txt", "w") as f1:
    f1.write(text)