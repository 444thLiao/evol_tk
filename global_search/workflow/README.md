
# workflow For dating

## 1. Download data

`ncbi-genome-download -s genbank -F fasta,genbank -A assembly_ids.list -p 20 -N bacteria`


## 2. postdownload (annotate the genomes or extract proteins sequences)
`python3 ~/script/evolution_relative/global_search/workflow/postdownload.py ./genbank ./genome_protein_files`


## 3. Run Orthofinder For building a genome/species tree

`orthofinder -f raw_genome_proteins -og -a 50 -1 -s diamond`

## 4. 