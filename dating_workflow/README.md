
# workflow For dating

## 1. Download data

`ncbi-genome-download -s genbank -F fasta,genbank -A assembly_ids.list -p 20 -N bacteria`


## 2. postdownload (annotate the genomes or extract proteins sequences)
`python3 /home-user/thliao/script/evolution_relative/global_search/workflow/postdownload.py ./genbank ./genome_protein_files`


## 3. Run Orthofinder For building a genome/species tree

`orthofinder -f raw_genome_proteins -og -a 50 -1 -s diamond`


## 4. extract conserved protein to following dating
`python3 /home-user/thliao/script/evolution_relative/ForOrthofinder/api/getSeqofOG_pro.py -i /home-backup/thliao/nitrification_for/dating_for/raw_genome_proteins/OrthoFinder/Results_Nov09_1/Orthogroups/Orthogroups.tsv -o ./og_extracted -rr 0.7 -doMSA`

`for f in `ls *.faa`; do mafft --maxiterate 1000 --genafpair --thread -1 > $.fa; done `
## 5. alignment these step4 

`python3 /home-user/thliao/script/evolution_relative/ForOrthofinder/api/concat_aln.py -i ./og_extracted_90`

## 6. build the tree and cluster into smaller tree
`iqtree -nt 50 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre og_extracted_90/iqtree -s og_extracted_90/concat_aln.aln -spp og_extracted_90/concat_aln.partition`

`/home-user/thliao/script/TreeCluster/TreeCluster.py -i  og_extracted_90/iqtree.treefile -tf argmax_clusters -t 1 > og_extracted_90/iqtree_treecluster`

`python3 /home-user/thliao/script/evolution_relative/dating_workflow/step_script/subset_tre_with_cluster.py -i  og_extracted_90/iqtree_treecluster -o og_extracted_90/new_genomes.list -s 100`

## 6. extract conserved proteins


## 7. convert it to phylip formatted file
`python3 /home-user/thliao/script/evolution_relative/dating_workflow/step_script/aln2phy.py -i "*.aln" `
``

## 8. Rough estimation of the substitution rate

