# Overview 

This repository is mainly designed for performing evolutional analysis quickly as you can.

Up to now, there are several parts for evolutional analysis.

Be aware of the following descriptions, the `$EVO` mean the base directory of this repository.

You could declare the environmental variable with a command like this `export EVO='$HOME/script/script_evolution' `

## 1. Toolkit
**The first and the most diversified bunches of scripts.**

In details, you could see specific document within `$EVO/bin` and `$EVO/api_tools`. (not done right now.)

Within the kit(`$EVO/bin`), there also have some different functional scripts. 
* batch script (for support some pipelines, it also contains a useful script like `batch_any.py` which could support some basically usage.) 
* convertor for NCBI (it could help to convert different id within the NCBI system. Like converting(**matching back**) protein id into the genome id. )
* convertor of others system (like *kegg* for now )
* download api from berry sequencing company (like `request_berry.py` with given token and directory need to be downloaded)
* script for **newick manipulation** `format_newick.py` (like reroot, add calibrations, sort etc...)

Within the kit(`$EVO/api_tools`), there also have some scripts. 

> The difference between `$EVO/bin` and `$EVO/api_tools` is the **usability**. Within `$EVO/bin`, it must a script easily to be used. But for `$EVO/api_tools`, normally it just some defined function which can only be used for pythoner. Of course, it have some exceptions due to the bad design at early period.

* iTOL formatting function (`itol_func.py` and `itol_template` which stodge modified template file download from [help page of iTOL](https://itol.embl.de/help.cgi))
* processing metadata from NCBI. Classifying the genomes into different habitats. (`metadata_for`)
* functions for formatting newick (`for_tree`)


This toolkit should continue updating during my working.

## 2. Post-analysis kit for specific Software

Within directory `for_software`, it will contain some software specific scripts. It should also used at other function depend on the validated IO.

Up to now, two kinds of software are well prepared for post-analysis like **badirate** and **bayestraits**

under `$EVO/for_software`

## 3. dating analysis specific scripts

These Scripts and some external scripts within toolkit could help ones to finish whole dating pipelines. 

For more details, see `doc/dating.md`

## 4. scripts for Orthofinder (also for pan-genome analysis)

These scripts for now is totally messy. It have two main script could easily be used.

* split out the orthologs which classified by orthofinder (it will be described clearly under doc `` )
* resort the orthologs after splitting with a reference genome.

## 5. scripts for aliTV (easily to see the pan-genome similarity visually ) 
incoming...