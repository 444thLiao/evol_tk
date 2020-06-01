# For processing Orthofinder output



## Part description
### api

This directory is the major part which contains script for processing Orthofinder output.

Simple introduction of these scripts:

1. getSeqofOG.py

    **Extracting genes from the output of Orthofinder**
    
    This provide multiple scenario for people to use it.
    
    * manually required some OG
    
        `python3 getSeqofOG.py -i Orthogroups.csv -og OG00001,OG00002 -o output_dir -single`
        
    * required all single copy and all presence OG
    
        `python3 getSeqofOG.py -i Orthogroups.csv -o output_dir -single -all_g`
        
    * required OG which is presenced at more than `100` genomes
    
        `python3 getSeqofOG.py -i Orthogroups.csv -o output_dir -single -t 100`
        
    * required OG which is presenced at specific genomes(with a input genome list file)
    
        `python3 getSeqofOG.py -i Orthogroups.csv -o output_dir -single -i2 genome_list.txt`
        
2. resort_OG_with_gaps.py & split_out_duplicated.py

> split out duplication within single Orthologs  
> 
> `-i`: Normally, the `Orthogroups.csv` with locus name
> 
> `-o`: output file
> 
> `-p`: normally, the prokka output directory. But considerate the flexible input. You could also pass a pattern like './\*/*.gbk' and parameter `-use_pattern`
> 
> `-gbk`: normally, it use gff file to retrieve genomic information. Pass this parameter to enable the utilization of genbank file
> 
> `-t`: threads to parse files. default is 20
> 
> `-use_pattern`: use pattern to retrieve genbank/gff files

    * split out the duplicated genes within one Orthogroups
        `python ~/script/evolution_relative/ForOrthofinder/bin/split_out_duplicated.py -i ./ortho_test/data/Results_Aug09/Orthogroups.csv -o ./test_splitted.csv -p ./prokka_o/ `

> resort_OG_with_gaps.py:
> 
> `-i`: Normally, the `Orthogroups.csv` with locus name after above script
> 
> `-o`: output sorted files
> 
> `-bc`:which columns you want to taken as backbone. default is the first one. Above script will sort the columns with the number of contigs. The genome with smaller number of contigs will ascend at the leftmost position.
> 
> `-sub`: Accept a file with genome names which is identical to the column in the table

* resort the orthlogroups table
`python ~/script/evolution_relative/ForOrthofinder/bin/resort_OG_with_gaps.py -i ./test_splitted.csv -o ./resorted.csv `


## raw

Some raw script for processing the data

## vis

Some visualization script.
For now, 
* it contains a script and some templates for iTOL
* a script drawing the genes summary bar plot
