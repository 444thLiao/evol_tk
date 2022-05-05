import io
import json
import os
from os.path import abspath
from os.path import exists, dirname, join, expanduser

import pandas as pd
from Bio import Entrez
from ete3 import NCBITaxa
from tqdm import tqdm

from api_tools.third_party.EntrezDownloader import EntrezDownloader
from api_tools.third_party.metadata_parser import (
    parse_bioproject_xml,
    parse_biosample_xml,
    parse_sra_xml,
    parse_assembly_xml,
)

import Bio


def eread(x):
    if float(Bio.__version__) >= 1.77:
        return Entrez.read(io.BytesIO(x.encode()))
    else:
        return Entrez.read(io.StringIO(x))


ncbi = NCBITaxa()
taxons = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
tmp_dir = "~/.tmp_getINFO"
tmp_dir = expanduser(tmp_dir)

import hashlib

email = os.getenv("EMAIL")
api_key = os.getenv("EKEY")
edl = EntrezDownloader(
    # An email address. You might get blocked by the NCBI without specifying one.
    email=email,
    # An API key. You can obtain one by creating an NCBI account. Speeds things up.
    api_key=api_key,
    num_threads=30,  # The number of parallel requests to make
    # The number of IDs to fetch per request
    batch_size=50,
    pbar=True,  # Enables a progress bar, requires tqdm package
)


def shash(astr):
    return hashlib.md5(astr.encode()).hexdigest()


def access_intermedia(obj, suffix="", redo=False):
    """
    It is mainly stodge intermediate result from previous results
    or provide a way to access it.
    Normally it should archive into current directory with md5 hashed file name.

    obj is necessary, normally it is a list of IDs/dictionary which also needed to generate md5 hashed file name.
    ofile is optional depend on what you want
    """
    nameofid = list(set(obj))
    _md5 = str(shash(";".join(list(sorted(nameofid)))))
    ofile = join(tmp_dir, _md5) + "_" + suffix
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile), exist_ok=1)
    if redo and exists(ofile):
        os.remove(ofile)
    if exists(ofile) and not redo:
        load_obj = json.load(open(ofile, "r"))
        if isinstance(load_obj, dict):
            tqdm.write("Dectect same cache, use it instead of run it again.")
            return load_obj
    else:
        if isinstance(obj, list):
            # just used to validate run it yet?
            pass
        elif isinstance(obj, dict):
            json.dump(obj, open(ofile, "w"))
        else:
            raise Exception("unexpected data type")


def parse_id(infile, columns=0):
    """
    1. format protein accession ID instead of using wrong one
    2. retrieving information provided inside the infile
    """
    id_list = []
    id2info = {}
    header_info = ""
    for row in tqdm(open(infile, "r")):
        if row.startswith("#"):
            # header row
            header_info = row.strip("\n").strip("#").split("\t")
            continue
        if row:
            id = row.split("\t")[columns].strip().strip("\n")
            if "||" in id:
                id = id.split("||")[-1]
            id_list.append(id)
            if not header_info:
                id2info[id] = ";".join(row.split("\t")[columns + 1 :]).strip("\n")
            else:
                id2info[id] = dict(
                    zip(header_info[1:], row.strip("\n").split("\t")[columns + 1 :])
                )
    return id_list, id2info


def unpack_gb(prot_t):
    """
    For standardizing unpack process from genbank formatted file
    """
    if not prot_t:
        return {}
    _cache_dict = {}
    # source feature (some have been implemented into annotations)
    source_fea = [_ for _ in prot_t.features if _.type == "source"]
    if source_fea:
        source_fea = source_fea[0]
        v = source_fea.qualifiers.get("isolation_source", [])
        v = [_ for _ in v if _]
        if v:
            _cache_dict[f"isolation_source"] = ";".join(v)
        v = source_fea.qualifiers.get("db_xref", [])
        v = [_ for _ in v if _ and _.startswith('taxon')]
        if v:
            _cache_dict[f"taxon"] = int(v[0].split(":")[-1])
    # above annotations
    annotations = prot_t.annotations
    ref_texts = [
        _
        for _ in annotations.get("references", [])
        if "Direct" not in _.title and _.title
    ]
    seq = str(prot_t.seq)
    org = annotations.get("organism", "")
    source = annotations.get("source", "")
    db_ = dict([_.split(":") for _ in prot_t.dbxrefs if ":" in _])
    biop = db_.get("BioProject", "")
    bios = db_.get("BioSample", "")
    nuccore_ID = annotations.get("db_source", "").split(" ")[-1]
    kw = ";".join(annotations.get("keywords", []))
    comment = annotations.get("comment", "")
    for idx in range(10):
        if idx < len(ref_texts):
            ref_t = ref_texts[idx]
            _cache_dict[f"reference title {idx + 1}"] = ref_t.title
            _cache_dict[f"reference journal {idx + 1}"] = ref_t.journal
            _cache_dict[f"reference authors {idx + 1}"] = ref_t.authors
    _cache_dict["organism"] = org
    _cache_dict["sequence"] = seq
    _cache_dict["source"] = source
    _cache_dict["bioproject"] = biop
    _cache_dict["biosample"] = bios
    _cache_dict["nuccore ID"] = nuccore_ID
    _cache_dict["keyword"] = kw
    _cache_dict["comment"] = comment
    return _cache_dict


def tax2tax_info(taxid):
    lineage = ncbi.get_lineage(int(taxid))
    rank = ncbi.get_rank(lineage)
    rank = {v: k for k, v in rank.items()}
    names = ncbi.get_taxid_translator(lineage)
    _dict = {}
    for c in taxons:
        if c in rank:
            _dict[c] = names[rank[c]]
    return _dict


def process_path(path):
    if not "/" in path:
        path = "./" + path
    if path.startswith("~"):
        path = expanduser(path)
    if path.startswith("."):
        path = abspath(path)
    return path


def parse_ipg(t):
    """
Id      Source  Nucleotide Accession    Start   Stop    Strand  Protein Protein Name    Organism        Strain  Assembly
13869719        RefSeq  NC_009937.1     1166617 1168176 -       WP_012169570.1  nitrogenase molybdenum-iron protein subunit beta        Azorhizobium caulinodans ORS 571        ORS 571      GCF_000010525.1
13869719        RefSeq  NZ_SOBM01000006.1       312912  314471  +       WP_012169570.1  nitrogenase molybdenum-iron protein subunit beta        Azorhizobium sp. AG788  AG788   GCF_004364705.1
13869719        INSDC   AP009384.1      1166617 1168176 -       BAF87037.1      nitrogenase molybdenum-iron protein beta chain  Azorhizobium caulinodans ORS 571        ORS 571 GCA_000010525.1
13869719        INSDC   SOBM01000006.1  312912  314471  +       TDT92758.1      Mo-nitrogenase MoFe protein subunit NifK        Azorhizobium sp. AG788  AG788   GCA_004364705.1
        :param t:
        :return:
    """
    headers = ['Id', 'Source', 'Nucleotide Accession', 'Start', 'Stop', 'Strand',
       'Protein', 'Protein Name', 'Organism', 'Strain', 'Assembly']
    bucket = []
    whole_df = pd.read_csv(io.StringIO(t), sep="\t", header=0)
    gb = whole_df.groupby('Protein')  # groupby the locus id
    all_dfs = [gb.get_group(x) for x in gb.groups]
    tqdm.write(f'processing {len(all_dfs)} dataframe retrieved')
    for indivi_df in all_dfs:
        indivi_df.index = range(indivi_df.shape[0])
        # aid = indivi_df.iloc[0, 6]
        indivi_df = indivi_df.fillna("")
        pos = []
        gb = []
        for _idx,row in indivi_df.iterrows():
            pos.append((row['Start'], row['Stop'], row['Strand'], row['Nucleotide Accession']))
            gb.append(row['Assembly'])
        for _idx,row in indivi_df.iterrows():
            bucket.append((row['Protein'], pos, gb))
    return bucket


def get_GI(ids, db, edl, no_order=False):
    tqdm.write("Get GI......")
    ids = [str(_) for _ in ids]
    if no_order:
        results, failed = edl.esearch(
            db=db, ids=ids, result_func=lambda x: eread(x)["IdList"], batch_size=50
        )
        all_GI = [_[-1] if type(_) == tuple else _ for _ in results]
        return all_GI
    results, failed = edl.esearch(
        db=db, ids=ids, result_func=lambda x: eread(x)["IdList"], batch_size=1
    )

    # for self.edl.esearch, it will auto **zip** searched term and its result.
    id2gi = dict(results)
    id2gi = {pid: id2gi.get(pid, "") for pid in ids}
    # extra
    if db == "assembly":
        remained_ids = [
            str(_).replace("GCF", "GCA")
            for _ in ids
            if str(_).startswith("GCF") and _ not in id2gi
        ]
        if remained_ids:
            results, failed = edl.esearch(
                db=db,
                ids=remained_ids,
                result_func=lambda x: eread(x)["IdList"],
                batch_size=1,
            )
            id2gi.update(dict(results))
            id2gi.update({k.replace("GCA", "GCF"): v for k, v in dict(results).items()})
    return id2gi, failed


def get_assembly(assembly_id=None, all_GI=None):
    if all_GI is None:
        all_GI = get_GI(assembly_id, "assembly", edl, no_order=True)

    results, failed = edl.esummary(
        db="assembly", ids=all_GI, result_func=lambda x: parse_assembly_xml(x)
    )
    dbsummary = {}
    for _dict in results:
        _dict = dict(_dict)
        for aid, info_dict in _dict.items():
            info_dict["TaxId"] = info_dict["SpeciesTaxid"]
            dbsummary.update(info_dict)
    return dbsummary


def get_srr(srr_list=None, all_GI=None):
    if all_GI is None:
        all_GI = get_GI(srr_list, "sra", edl, no_order=True)
    # if len(srr_list) != len(all_GI):
    #     print(f"{len(results)-len(srr_list)} item might failed to retrieve its GI")

    results, failed = edl.efetch(
        db="sra",
        ids=all_GI,
        retmode="xml",
        retype="xml",
        batch_size=50,
        result_func=lambda x: parse_sra_xml(x),
    )
    srr2info = {}
    for _ in results:
        if isinstance(_, dict):
            srr2info.update(_)
    return srr2info


def get_bioproject(bp_list=None, all_GI=None):
    if all_GI is None:
        all_GI = get_GI(bp_list, "bioproject", edl, no_order=True)
    results, failed = edl.efetch(
        db="bioproject",
        ids=all_GI,
        retmode="xml",
        retype="xml",
        result_func=lambda x: parse_bioproject_xml(x),
    )
    bp2info = {}
    for _ in results:
        if isinstance(_, dict):
            bp2info.update(_)
    return bp2info


def get_biosample(bs_list=None, all_GI=None):
    if all_GI is None:
        all_GI = get_GI(bs_list, "biosample", edl, no_order=True)
    results, failed = edl.efetch(
        db="biosample",
        ids=all_GI,
        retmode="xml",
        retype="xml",
        result_func=lambda x: parse_biosample_xml(x),
    )
    bs2info = {}
    for _ in results:
        if isinstance(_, dict):
            bs2info.update(_)
    return bs2info


from Bio import SeqIO


def get_seq_with_pos(nid2pos):
    # nid2pos = {}
    # for idx,row in tqdm(sub_df.iterrows()):
    #     nid = row['nucleotide id']
    #     v = (row['start'],row['end'],row['strand'])
    #     if nid not in nid2pos:
    #         nid2pos[nid] = [v]
    #     else:
    #         nid2pos[nid].append(v)
    get_records = []
    for nid, v_d in tqdm(nid2pos.items()):
        for s, e, strand in v_d:
            s, e = int(s), int(e)
            if str(s) == "nan":
                continue
            strand = 1 if strand == "+" else 2
            results, failed = edl.efetch(
                db="nucleotide",
                ids=[nid],
                result_func=lambda x: list(SeqIO.parse(io.StringIO(x), "fasta")),
                batch_size=1,
                retype="fasta",
                retmode="text",
                seq_start=s,
                seq_stop=e,
                strand=strand,
            )
            get_records.extend(results)
    return get_records


import requests
from bs4 import BeautifulSoup


def hard_retrieval():
    # for some protein ids, it can not retrieve via Entrez but can access with web
    requests.get(f"https://www.ncbi.nlm.nih.gov/protein/pids")

    record = requests.get("https://www.ncbi.nlm.nih.gov/protein/MBI3798512.1")
    soup = BeautifulSoup(record.content)
    pass
