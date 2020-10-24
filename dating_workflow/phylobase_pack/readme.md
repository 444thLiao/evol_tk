# instruction for usage of phylobayes

> following instructions are suitable for phylobayes 4.1c

## practice example (here use plancto)
> mkdir run2; pb -d ../dating_for/phy_files/83g_concat.phy -T ./tree.newick -cal cal_set14.txt -lg -ugam -dgam 4 -sb 0.05 -bd -catfix C10 run2/test &
> mkdir run1; pb -d ../dating_for/phy_files/83g_concat.phy -T ./tree.newick -cal cal_set14.txt -lg -ugam -dgam 4 -sb 0.05 -bd -cat -catfix C10 run1/test &
> mkdir prior_only; pb -d ../dating_for/phy_files/83g_concat.phy -T ./tree.newick -cal cal_set14.txt -lg -ugam -sb 0.05 -bd -prior -ncat 1 -dgam 1 prior_only/test &
> readdiv -x 10000 run1/test & ; readdiv -x 10000 run2/test & ; readdiv -x 10000 prior_only/test &

## Input data format
Aligned Sequences:
    PHYLIP format

Trees:
    NEWICK format

Outgroups

Calibrations (Upper or lower limits can be set equal to -1, in which case no limit is enforced.)

## running a chain
pb -d <dataset> <chainname>


## stop a chain
### soft-stop a chain
stoppb <chainname>  

or 

echo 0 > <chainname>.run

## restart a chain
pb <chainname>


## Checking convergence
1. visualize trace file(<chainname>.trace) via **Tracer**
2. bpcomp -x 1000 10 <run1> <run2>   (burn-in of 1000, and sub-sampling every 10 trees,)
Some guidelines:
   * maxdiff < 0.1: good run.
   * maxdiff < 0.3: acceptable: gives a good qualitative picture of the posterior consensus.
   * 0.3 < maxdiff < 1: the sample is not yet suffciently large, and the chains have not
   converged, but this is on the right track.
   * if maxdiff = 1 even after 10,000 points, this indicates that at least one of the runs is
   stuck in a local maximum.
3. tracecomp -x 1000 <chain1> <chain2>
   * maxdiff < 0.1 and minimum effective size > 300: good run;
   * maxdiff < 0.3 and minimum effective size > 50: acceptable run.


## Obtaining posterior consensus trees and parameter estimates
readpb -x 1000 10 <chainname>
    * <chainname>_sample.general, containing a few general statistical averages (mean
    posterior log likelihood, tree length, alpha parameter, number of modes, etc.)
    * <chainname>_sample.con.tre : the majority-rule posterior consensus tree
    * <chainname>_sample.bplist : the list of weighted bipartitions.
  
## Divergence times
### read the div times
readdiv -x 10000 <chainname>

## parameters would be used (for pb)
-ln: log-normal autocorrelated relaxed clock
-rp: prior time of the root
-sb: soft bounds/constraints (which works only under the birth death prior)














