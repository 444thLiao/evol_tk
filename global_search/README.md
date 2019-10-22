# general guidelines

## Purpose:
For retrieving protein sequence relative to interested genes.
1. retrieve sequence
2. retrieve relative infomation from biosample and bioproject


## pipelines
1. `diamond` search from `nr` database according to `prepare_db.py`
2. parse diamond output into a long list by `parse_diamond_result.py`
3. parse again to get sequence from `nr` by `simple_get_fa_from_diamond.py`
4. using `TIGFAM` and `kegg` to double check if this sequence is requested gene by `hmmsearch` according to 
5. use `bin/fasta2id_list.py` to convert filtered fasta file into a idlist file.
6. use `get_info.py` to retrieve relative information from NCBI by **Entrez**.
