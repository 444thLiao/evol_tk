# For processing Orthofinder output



## Part description
### api

This directory is the major part which contains these script for processing Orthofinder output.

Below I will simple introduce these scripts

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

## raw

Some raw script for processing the data

## vis

Some visualization script.
For now, 
* it contains a script and some templates for iTOL
* a script drawing the genes summary bar plot
