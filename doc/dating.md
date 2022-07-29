# dating analysis with written scripts
The doc is introducing a pipelines from downloaded genomes to molecular dating analysis with embedded scripts.


In beginning, declare a **env variable** for usage.

> `export EOV=/home-user/thliao/script/evolution_relative`

## 1. Download data
check the taxonomy id from https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi

Like taxid of phylum Cyanobacteria is **1117**, Of course you could use some smaller lineage like Bradyrhizobium **374**.

`export NAME=cyano; export NCBIDB=~/.cache/ncbi-genome-download/refseq_bacteria_assembly_summary.txt;`

`python3 ~/software/ncbi-genome-download/contrib/gimme_taxa.py -o ./$NAME_phylum.txt 1117;`

`cat ./$NAME_phylum.txt | cut -f2 > only_tid.txt;`

`cut -f1,6,7 $NCBIDB | grep -wFf only_tid.txt | cut -f1 > ./assembly_ids.list;`

`ncbi-genome-download -s refseq -F fasta,genbank -A assembly_ids.list -p 20 -N bacteria;`

> **NOTE**
> 
> The file used for getting the genome id should be download from [genbank](ftp://ftp.ncbi.nih.gov/genomes/genbank/assembly_summary_genbank.txt) or [refseq](ftp://ftp.ncbi.nih.gov/genomes/refseq/assembly_summary_refseq.txt)
> 
> In above example, it used **~/.cache/ncbi-genome-download/refseq_bacteria_assembly_summary.txt**. So at the following command, you should pass **refseq** to `-s`. Or it will get wrong.


## 2. postdownload (annotate the genomes or extract proteins sequences)
`python3 $EOV/dating_workflow/toolkit/postdownload.py ./genbank ./genome_protein_files`

If you want to strictly follow above commands, you could also pass nothing to this script.

positional arguments : first is the downloaded directory (it  )
second is the output dir (it will auto stodge the proteins sequence into it and **format** the file name)

> **Formatting Rule**
>
> if you pass GCA_123456.2, it will be converted into A123456v2. If you pass GCF_123456.1, it will be converted into F123456v1.
>
> With this simply conversion, it could also easily be converted back. 


Due to the differences of name conversion between genbank and refseq, it need to pass -f

## 3.1 Run Orthofinder For building a genome/species tree

For small number of genomes, you could use the result of orthofinder to build the phylogenetic tree.

`orthofinder -f raw_genome_proteins -og -a 50 -1 -S diamond`

## 3.2 extract conserved protein for building phylogenetic tree
### a. extract proteins
For larger number of genomes, running orthofinder is time wasted. It could use some reported conserved proteins to build phylogenetic tree. 

`python3 $EOV/dating_workflow/step_script/extract_bac120.py -in_p '../rawdata/genome_protein_files' -in_a ./bac120_annotate -o ./bac120_extract/seq_e30 -evalue 1e-30`

`python3 $EOV/dating_workflow/step_script/extract_cog25.py -in_p '../rawdata/genome_protein_files' -in_a ./cog25_annotate -o ./cog25_single/seq_e20 -evalue 1e-20`

`python3 $EOV/dating_workflow/step_script/extract_r50.py -in_p '../rawdata/genome_protein_files' -in_a ./r50_annotate -o ./r50_single/seq_e20 -evalue 1e-20`

* in_p: input proteins directory
* in_a: actually an output annotation directory
* o: output directory contains series of proteins fasta.
* evalue: cutoff evalue for extracting validated proteins after annotating (doesn't use it during annotating.)

Up to day, there are three sets of proteins/genes could be taken as target proteins for following phylogenetic analysis.
The parameters of them are similar.

### b. alignment, trimal and concat them into seq with partition.
Following commands using 187 genomes to perform analysis

`python3 ~/bin/batch_run/batch_mafft.py -i ./cog25_single/seq_e20 -s faa -o ./cog25_single/187g_aln -f -m ginsi -gl ./dating_for_187g.list;`

`python3 ~/bin/batch_run/batch_trimal.py -i ./cog25_single/187g_aln -o ./cog25_single/187g_aln -ro 0.5 -so 40;`  

> Adjusting the `ro` and `so` for ensuring the existing of backbone gene.
>
> ro: **-resoverlap**,Minimum overlap of a positions with other positions in the column to be considered a "good position". Range: [0 - 1]. (see User Guide).
>
> so: **-seqoverlap**,Minimum percentage of "good positions" that a sequence must have in order to be conserved. Range: [0 - 100](see User Guide).

`python3 ~/script/evolution_relative/dating_workflow/toolkit/concat_aln.py -i ./cog25_single/187g_aln -ct phy -gl ./dating_for_195g.list -o ./dating_for/phy_files/187g_concat.trimal -s trimal -no_graph;`

### c. build tree
`iqtree -nt 50 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre og_extracted_90/iqtree -s og_extracted_90/concat_aln.aln -spp og_extracted_90/concat_aln.partition`

## 4. preparing tree for dating analysis

concating two tree into one. And most important one is adding calibrations. 

`format_newick.py cat -i ~/data/cyano_basal/analysis/pruned_formatted.newick -i2 ./trees/final/195g_cyano_complex_bac120.formatted.newick -o ./trees/final/198g_merged.newick -f 3 -f_to 3 `

`format_newick.py add-cal -i ./trees/final/198g_merged.newick -c ./dating_for/calibration2.txt ./dating_for/198g_3cal.newick  `

## 5. performing dating analysis
`python3 ~/script/evolution_relative/dating_workflow/step_script/dating_pro.py -i ./dating_for/phy_files/83g_concat.phy -it ./dating_for/cal_tree/83g_set1.newick -o ./dating_for/83g/83g_set1/ -p 1 -rg '1 100 1' -sf 10 -c 3`

details of parameters could use `python3 script.py --help` to get help.


## 6. process the result into iTOL suitable file
`python3 ~/script/evolution_relative/dating_workflow/figtree2itol.py -i ./trees/iqtree/168g_concat.treefile -i2 ./dating_for/168g_bac120_3cal/mcmc_for/FigTree.tre -o ./dating_for/result_tree/168g_bac120_3cal_dating.newick -r 'GCA_000011385.1,GCA_000332175.1' `

The first tree could get bootstrap values for annotation. The second one is mainly annotated with time.
