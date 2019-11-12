


seqfile = 00_seqs.phy
treefile = 03_iqtree.point.tre
outfile = 03_mlc.txt

ndata = 1
seqtype = 2     * amino acid sequence

noisy = 3       * how much rubbish on the screen
verbose = 0     * concise output
runmode = 0     * user tree

model = 2       * model of amino acid substitution (empirical model, given in the file aaRatefile)
aaRatefile = lg.dat     *empirical amino acid substitution rate matrix
Mgene = 0

fix_alpha = 0
alpha = 0.5
Malpha = 0
ncatG = 5
fix_rho = 1
rho = 0.

clock = 1       * strick clock, see PAML reference
getSE = 1
cleandata = 0
RateAncestor = 0        *
