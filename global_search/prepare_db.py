
# latest database including nr,env_nr,swissprot



# run diamond 
query_fa = "all.faa"
nr_db = "/home-user/thliao/data/protein_db/NCBI/nr"
cmd_template = f"diamond blastp -q {query_fa} -d  -o query_result/nr_retrieve_all.diamond -k 0 -e 1e-20 --max-hsps 1 -f 6 qseqid sseqid salltitles pident length evalue bitscore"
# params explanation
# evalue choose 1e-20 to filter out false positive 
# k choose 0 to get all aligned result instead of default 25
# max-hsps 1 for just remained the highest one
# output format 6 for formatting the out result 