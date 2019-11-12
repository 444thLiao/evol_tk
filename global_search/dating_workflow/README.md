
# workflow For dating

## 1. Download data

`ncbi-genome-download -s genbank -F fasta,genbank -A assembly_ids.list -p 20 -N bacteria`


## 2. postdownload (annotate the genomes or extract proteins sequences)
`python3 ~/script/evolution_relative/global_search/workflow/postdownload.py ./genbank ./genome_protein_files`


## 3. Run Orthofinder For building a genome/species tree

`orthofinder -f raw_genome_proteins -og -a 50 -1 -s diamond`


## 4. extract conserved protein to following dating
`python3 ~/script/evolution_relative/ForOrthofinder/api/getSeqofOG_pro.py -i /home-backup/thliao/nitrification_for/dating_for/raw_genome_proteins/OrthoFinder/Results_Nov09_1/Orthogroups/Orthogroups.tsv -o ./og_extracted -rr 0.7 -doMSA`


## 5. alignment these step4 

`python3 ~/script/evolution_relative/ForOrthofinder/api/concat_aln.py -i ./og_extracted_90`

## 6. extract conserved proteins
```

```
