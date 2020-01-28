# comands of pipelines instructions in brief

## download a phylum of genome data from NCBI
check the taxonomy id from 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi'
etc.  taxid of phylum Cyanobacteria is 1117

`python3 ~/software/ncbi-genome-download/contrib/gimme_taxa.py -o ./cyano_phylum.txt 1117`

`cat ./cyano_phylum.txt|cut -f2 > only_tid.txt`

`cut -f1,6,7 ~/.cache/ncbi-genome-download/refseq_bacteria_assembly_summary.txt | grep -wFf only_tid.txt | cut -f1 > ./assembly_ids.list`

`ncbi-genome-download -s refseq -F protein-fasta -A assembly_ids.list -p 20 bacteria`

if you want to download from genbank(more), switch **refseq into genbank** 
## if you want to retrieve metadata of these genomes, you could...
`python3 ~/script/evolution_relative/bin/ncbi_convertor/pid2bio.py -i ./assembly_ids.list -o ./biometadata.csv -s genome`

## pack up all proteins
`python3 ~/script/evolution_relative/dating_workflow/toolkit/postdownload.py ./refseq ./genome_protein_files`
positional arguments : first is the downloaded dir (dependent on the name after `-s` at previous)
second is the output dir (it will auto stodge the proteins sequence into it and format the file name)

## extract target proteins for following phylogenetic analysis (either)
`python3 /home-user/thliao/script/evolution_relative/dating_workflow/step_script/extract_bac120.py -in_p '../rawdata/genome_protein_files' -in_a ./bac120_annotate -o ./bac120_extract/seq_e30 -evalue 1e-30`

`python3 /home-user/thliao/script/evolution_relative/dating_workflow/step_script/extract_cog25.py -in_p '../rawdata/genome_protein_files' -in_a ./cog25_annotate -o ./cog25_single/seq_e20 -evalue 1e-20`

`python3 /home-user/thliao/script/evolution_relative/dating_workflow/step_script/extract_r50.py -in_p '../rawdata/genome_protein_files' -in_a ./cog25_annotate -o ./cog25_single/seq_e20 -evalue 1e-20`

Up to day, there are three sets of proteins/genes could be taken as target proteins for following phylogenetic analysis.

The parameters of them are similar.

## construct tree with above target proteins
export extract_dir

`python3 ~/bin/batch_run/batch_mafft.py -i ./bac120_extract/seq_e50 -s faa -o ./bac120_extract/232g_cyano -f -m einsi -gl ../rawdata/assembly_ids.list -fix_ref`
`python3 ~/bin/batch_run/batch_trimal.py -o ./bac120_extract/232g_cyano -i ./bac120_extract/232g_cyano`
`python3 /home-user/thliao/script/evolution_relative/dating_workflow/toolkit/concat_aln.py -i ./bac120_extract/232g_cyano -o ./trees/concat/cyano_concat.trimal -s trimal -gl ../rawdata/assembly_ids.list -ct partition -no_graph`
`iqtree -nt 35 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre ./trees/iqtree/cyano_concat -s ./trees/concat/cyano_concat.trimal -spp ./trees/concat/cyano_concat.partition`
`iqtree -nt 30 -wbtl -bb 1000 -m LG+C20+F+G -redo -ft ./trees/iqtree/cyano_concat.contree -s ./trees/concat/cyano_concat.trimal -pre ./trees/iqtree/cyano_complex_bac120`
