# comands of pipelines instructions in brief

## download a phylum of genome data from NCBI
check the taxonomy id from 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi'
etc.  taxid of phylum Cyanobacteria is 1117
`python3 ~/software/ncbi-genome-download/contrib/gimme_taxa.py -o ./cyano_phylum.txt 1117`

`cat ./cyano_phylum.txt|cut -f2 > only_tid.txt`

`cut -f1,7 ~/.cache/ncbi-genome-download/refseq_bacteria_assembly_summary.txt | grep -wFf only_tid.txt | cut -f1 > ./assembly_ids.list`

`ncbi-genome-download -s refseq -F fasta -A assembly_ids.list -p 20 bacteria`
if you want to download from genbank(more), switch **refseq into genbank** 

## further analysis
`python3 ~/script/evolution_relative/dating_workflow/toolkit/postdownload.py ./refseq ./genome_protein_files`