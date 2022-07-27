from os.path import *
from tqdm import tqdm
from glob import glob
import os
import pandas as pd
from ete3 import Tree
from subprocess import check_call
from collections import defaultdict,Counter
HMM_DBPATH = "/mnt/home-backup/thliao/kofam/kegg.HMM"
KEGG_AnnoODIR = "/mnt/maple/thliao/data/NCBI/modified_data/annotations/kegg"
NCBI_PROTEIN = "/mnt/maple/thliao/data/NCBI/modified_data/direct_protein_files"
NCBI_PROKKA = "/mnt/maple/thliao/data/NCBI/modified_data/direct_protein_files"

def kegg_anno_cmds(aids,cpu=32):
    no_faa_gids = []
    # ! kegg annotation
    kegg_all_hmm = HMM_DBPATH
    target_dir = KEGG_AnnoODIR
    cmds = []
    for gid in tqdm(aids):
        faa = f"{NCBI_PROTEIN}/{gid}.faa"
        prokka_faa = f"{NCBI_PROKKA}/{gid}/{gid}.faa"
        if not exists(faa) and not exists(prokka_faa):
            no_faa_gids.append(gid)
        elif not exists(faa) and exists(prokka_faa):
            os.system(f"ln -s {prokka_faa} {faa}")
        ofile = join(target_dir, f"{gid}.tab")
        cmd = f"hmmsearch --tblout {ofile} --acc --noali --notextw --cpu {cpu} {kegg_all_hmm} {faa} "
        if not exists(ofile):
            cmds.append(cmd)
    return cmds,no_faa_gids

def get_no_fna(aids):
    missing_aids = []
    curr_dir = join(db_dir, 'genbank','bacteria')
    sub_aids = [_ 
                for _ in tqdm(aids)
                if not glob(join(curr_dir, _, '*.fna.gz'))]
    missing_aids.extend(sub_aids)
    return missing_aids



def convert_genome_ID(genome_ID):
    # for GCA_900078535.2
    # it will return 900078535v2
    if isinstance(genome_ID, str) and genome_ID.startswith('GC'):
        return genome_ID.split('_')[-1].replace('.', 'v')
    else:
        return genome_ID


def convert_genome_ID_rev(locus_ID, prefix='GCA_',not_add_prefix_ids=[]):
    # for 900078535v2
    # it will return prefix + 900078535.2
    if locus_ID in not_add_prefix_ids:
        return locus_ID
    if '|' in str(locus_ID):
        # other labmate used
        genome_name = locus_ID.partition('|')[0]
        return genome_name

    if isinstance(locus_ID, str) and not locus_ID.startswith('GC'):
        if '_' in locus_ID:
            # tianhua version, it won't contain |
            locus_ID = locus_ID.partition('_')[0]
            if locus_ID in not_add_prefix_ids:
                return locus_ID
            else:
                return prefix + locus_ID.replace('v', '.')
        else:
            return prefix + locus_ID.replace('v', '.')
    else:
        return locus_ID


def run(cmd):
    check_call(cmd,
               shell=True,
               stdout=open('/dev/null', 'w'),
               )
def refresh_tmp(tmpdir):
    if not exists(tmpdir):
        os.system(f"mkdir -p {tmpdir}")
    else:
        os.system(f"rm -rf {tmpdir}/*")

def get_files(in_path,suffix):
    suffix = suffix.strip('.')
    files = []
    if ',' in in_path:
        in_path = in_path.split(',')
    else:
        in_path = [in_path]
    for _ in in_path:
        _f = glob(join(_, '*.' + suffix.strip('.')))
        files.extend(_f)
    return files
    
def get_genomes(genome_list,
                simple_concat=True):

    """
    Accepting a file. 
    It could only contain a column of genome names for simple jobs.
    
    Or
    It could contains multiple lines separated with TAB.
    Besides the first column, the following columns should be the gene names in the alignment.

    Returns:
        dict: name2grouping, maybe empty grouping
    """
    # if genome_list is None:
    #     genome_list = join(indir, 'selected_genomes.txt')
    if genome_list is None:
        return
    if ',' in genome_list:
        list_vals = [_.strip() for _ in genome_list.split(',')]
        if os.path.isfile(list_vals[0]):
            rows = []
            for l in list_vals:
                rows.extend(open(l, 'r').read().split('\n'))
        else:
            return {k:k for k in list_vals}  # simple id list
    else:
        rows = [_ for _ in open(genome_list, 'r').read().split('\n') if _ ]
    
    final_name2grouping = defaultdict(set)
    for row in rows:
        if '\t' not in row:
            name = row if simple_concat else convert_genome_ID(row)
            final_name2grouping[name].add(name)
        else:
            name = row.split('\t')[0]
            final_name2grouping[name].add(row.split('\t')[1])
    final_name2grouping = {k:v for k,v in final_name2grouping.items() if k}
    return final_name2grouping


def get_tophit(a2list_b, top_hit):
    a2list_b = {k: tuple(sorted(v,key=lambda x:x[1]))
                 for k,v in a2list_b.items()}  # ascending, nearly unchanged
    if top_hit:
        a2single_b = {k: [v[0][0]]
                     if v else []
                     for k, v in a2list_b.items()}
        b2_a = {}
        for a,list_b in a2list_b.items():
            for b,eval in list_b:
                if b in b2_a:
                    b2_a[b] = (a,eval) if b2_a[b][1]> eval else b2_a[b]
                else:
                    b2_a[b] = (a,eval)
    else:
        a2single_b = {k: [_[0] for _ in v]
                     if v else []
                     for k, v in a2list_b.items()}
        
        b2_a = defaultdict(list)
        for a,list_b in a2list_b.items():
            for b,eval in list_b:
                b2_a[b].append((a,eval))
        b2_a = {k: tuple(sorted(v,key=lambda x:x[1]))
                 for k,v in b2_a.items()}  # ascending, nearly unchanged
    return a2single_b,b2_a


def parse_blastp(ofile, match_ids=[], filter_evalue=1e-3, 
                 seq2length=None,
                 pos_location = [8,9,1,1,0], cov=80, 
                 top_hit=False):
    # pos_location should contain the index of specific cols
    # start, end, id of `seq2length`, unique id, mapped id.
    start,end,_pid,unique_id,mapped_id = pos_location
    if not match_ids:
        gid2locus = defaultdict(list)
    else:
        gid2locus = {k: [] for k in match_ids}
    for row in open(ofile, 'r'):
        sep_v = row.split('\t')
        locus = sep_v[mapped_id]
        evalue = float(sep_v[10])
        if seq2length:
            s,e = sep_v[start],sep_v[end]
            pid = sep_v[_pid]
            align_cov = abs(int(e)-int(s))/seq2length[pid]*100
            if align_cov <= cov:
                continue
        if filter_evalue and evalue > filter_evalue:
            continue
        if sep_v[unique_id] in match_ids:
            gid2locus[sep_v[unique_id]].append((locus, evalue))
        if not match_ids:
            gid2locus[sep_v[unique_id]].append((locus, evalue))
        # normally, the unique id should be the locus in the genome file. Thus, one locus must be only assgined to a gene (no need to adjust)
        # However, across all the locus list, we also need to select the best aligned locus to a gene (need to adjust)
        # Thus, the `top_hit` should be used to control top hit or not to the latter one(mapped id)
    a2single_b,b2_a = get_tophit(gid2locus, top_hit=top_hit)
    return a2single_b,b2_a

def parse_hmmscan(ofile, filter_evalue=1e-20, top_hit=False, _pos=[0,2],
                  evalue_pos=4):
    # _pos indicate the index of query and subject
    q_pos,s_pos = _pos
    gid2locus = defaultdict(list)

    for row in open(ofile, 'r'):
        if row.startswith('#'):
            continue
        r = row.split(' ')
        r = [_ for _ in r if _]
        gene_id = r[q_pos]
        locus_tag = r[s_pos]
        evalue = float(r[evalue_pos])
        if filter_evalue and evalue > filter_evalue:
            continue
        gid2locus[gene_id].append((locus_tag, evalue))
    a2single_b,b2_a = get_tophit(gid2locus, top_hit=top_hit)
    return a2single_b,b2_a




import io
def read_hmmsearch(domtblout,apply_preset_filter=True):
    # modified from https://github.com/aweimann/traitar
    # particuarly for domtblout
    # only apply for the output generated by the cmd f"hmmsearch --cut_ga --domtblout {ofile} --cpu 32 {kegg_all_hmm} {faa} " 
    # hmmer (v3.3)
    hmmer_colnames = ['target name','target accession','tlen','query name','accession','qlen','E-value','score_overall','bias_overall','#','of','c-Evalue','i-Evalue','score_domain','bias_domain','from_hmm','to_hmm','ali_from','ali_to','env_from','env_to','acc','description of target']
    cleaned = "".join(filter(lambda x: not x.startswith("#"), 
                             ["\t".join(i.split(None, 22)) 
                              for i in open(domtblout).readlines()]))
    m = pd.read_csv(io.StringIO(cleaned), sep = "\t",  header = None)
    m.columns = hmmer_colnames
    def filter_dbcan(m):
        """apply thresholds as suggested by the dbCAN authors i.e. a bit score threshold of 25. If the alignment is less than 80 bps long use an e-value threshold of 0.001 and or an evalue threshold of 0.00001 if the alignment is longer than 80 bps """
        return ((m.loc[:,"score_domain"] >= 25) & ((m.loc[:,"i-Evalue"] <= 0.001) & (m.loc[:,"ali_to"] - m.loc[:, "ali_from"] <= 80) | (m.loc[:,"i-Evalue"] <= 0.00001) & (m.loc[:,"ali_to"] - m.loc[:, "ali_from"] > 80))) 
    if apply_preset_filter:
        keep = filter_dbcan(m)
        final_df = m.loc[keep,:]
    else:
        final_df
    return final_df

def read_hmmsearch_tbl(tblout,apply_preset_filter=False):
    # modified from https://github.com/aweimann/traitar
    # particuarly for tblout
    # implement for "hmmsearch --tblout {ofile} --acc --noali --notextw --cpu 32 {kegg_all_hmm} {faa} "
    # hmmer (v3.3)
    hmmer_colnames = ['target name','target accession','query name','query accession',
                      'E-value(full)','score(full)','bias(full)','E-value(best1domain)','score(best1domain)','bias(best1domain)',
                      "exp","reg","clu","ov","env","dom","rep","inc","description of target"]
    cleaned = "".join(filter(lambda x: not x.startswith("#"), 
                             ["\t".join(i.split(None, 18)) 
                              for i in open(tblout).readlines()]))
    m = pd.read_csv(io.StringIO(cleaned), sep = "\t",  header = None)
    m.columns = hmmer_colnames
    def filter_dbcan(m):
        """apply thresholds as suggested by the dbCAN authors i.e. a bit score threshold of 25. If the alignment is less than 80 bps long use an e-value threshold of 0.001 and or an evalue threshold of 0.00001 if the alignment is longer than 80 bps """
        return ((m.loc[:,"score(full)"] >= 25) & (m.loc[:,"E-value(full)"] <= 1e-20)  ) 
    if apply_preset_filter:
        keep = filter_dbcan(m)
        final_df = m.loc[keep,:]
    else:
        final_df
    return final_df
    