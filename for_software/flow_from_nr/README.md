# general guidelines

## Purpose:
For retrieving protein sequence relative to interested genes.
1. retrieve sequence
2. retrieve relative infomation from biosample and bioproject


## pipelines
1. `diamond` search from `nr` database according to `prepare_db.py`
2. parse diamond output into a long list by `parse_diamond_result.py`
3. parse again to retrieve particular sequence from `nr` by `simple_get_fa_from_diamond.py`
4. using `TIGFAM` and `kegg` to double check if this sequence is requested gene by `hmmsearch` according to `filter_fa.py`
5. For gene have over 10K+ sequences, perform `cd-hit` first and get one seq from each cluster.
6. manually choose a suitable number sequence to go on.
6. use `bin/fasta2id_list.py` to convert filtered fasta file into a idlist file.
7. use `get_info.py` to retrieve relative information from NCBI by **Entrez**.

8. `mafft` and **build** tree
9. annotate tree result with `depostied` and `intermediated` dumpped file.