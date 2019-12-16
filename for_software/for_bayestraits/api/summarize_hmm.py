from tqdm import tqdm
from glob import glob
from collections import defaultdict
import pandas as pd
from dating_workflow.step_script import convert_genome_ID_rev
import click
from os.path import join


def retrieve_info(indir, suffix):
    gid2locus2ko = defaultdict(list)
    files_list = glob(join(indir, f'.{suffix}'))
    if not files_list:
        exit(
            f"no files could be found with input {join(indir,f'.{suffix}')},please check the parameters. ")
    tqdm.write("reading all annotated result")
    for hf in tqdm(files_list):
        for row in open(hf):
            if row.startswith('#'):
                continue
            r = row.split(' ')
            r = [_ for _ in r if _]
            gene_id = r[0]
            ko = r[2]
            evalue = float(r[4])
            gid2locus2ko[convert_genome_ID_rev(gene_id)].append(
                (gene_id, ko, evalue))
    return gid2locus2ko


def filtration_part(gid2locus2ko,evalue=1e-50):
    # filter out with hard threshold of evalue
    gid2locus2ko_e50 = {k: [(_[1], _[0], _[2])
                            for _ in v if _[2] <= evalue]
                        for k, v in tqdm(gid2locus2ko.items())}
    # select minimum evalue among all matched KO for each locus
    # TODO: it may be corrected at following
    used_locus = {}
    locus2ko = {}
    for gid, inlist in tqdm(gid2locus2ko_e50.items()):
        for key, v, evalue in inlist:
            if evalue <= used_locus.get(v, 100):
                used_locus[v] = evalue
                locus2ko[v] = key

    gid2locus2ko_e50 = defaultdict(lambda: defaultdict(list))
    for locus, ko in tqdm(locus2ko.items()):
        gid = convert_genome_ID_rev(locus)
        gid2locus2ko_e50[gid][ko].append(locus)

    gid2locus2ko_e50 = {g: {ko: ','.join(v) for ko, v in d.items()}
                        for g, d in gid2locus2ko_e50.items()}
    return gid2locus2ko_e50


tmp_df = pd.DataFrame.from_dict(gid2locus2ko_e50, orient='index')
tmp_df.to_csv(
    './protein_annotations/kegg_hmm_merged_info_e20.tab', index=1, sep='\t')


@click.command()
@click.option("-i", "indir",)
@click.option("-s", 'suffix')
@click.option("-e", "evalue")
@click.option("-t", "transpose")
def main(indir,suffix,evalue,transpose):
    gid2locus2ko = retrieve_info(indir,suffix)
    


if __name__ == "__main__":
    main()
