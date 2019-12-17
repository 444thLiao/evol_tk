# introduction of scripts under current directory

## API

For running bayestraits analysis with multistats+MCMC with input **tree** and **metadata**
`python3 ~/script/evolution_relative/for_software/for_bayestraits/api/run_habitat.py -i intree -o odir -im inmetadata`


For summarizing the results output from batch hmm(all versus each KO)

It is suitble for hmm output which query is KO fam, and the subject is your sequence.
`python3 ~/script/evolution_relative/for_software/for_bayestraits/api/summarize_hmm.py -i ./kegg_hmmsearch -o hmmsearch_merged -s tab`
