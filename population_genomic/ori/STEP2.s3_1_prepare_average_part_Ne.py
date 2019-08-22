#!/usr/bin/env python

import os
import sys

# We can now estimate the effective population size using neaverage.pl, which weighs samples by their recombination 
# distances, for which we need to inform it about the recombination files used to generate the file.

script = "STEP2.s3_2_average_part_Ne.sh"
ddir = "parts_Ne"
subset_interval = 5  # interval of subset size when perform EM inference of Ne

# get list of rfiles whose Ne was estimated in previous step
rfiles = []
for ffile in os.listdir(ddir):
    if ffile.endswith(".phase") and ffile.startswith("bialle_SNP.ref_"):
	rfiles.append(ffile.replace(".phase", ".recfile"))

# get number of individuals from one phase file
num_indv = 0
with open("%s/%s" % (ddir, rfiles[0].replace(".recfile", ".phase"))) as f:
    f.readline()
    num_indv = int(f.readline()[:-1]) 

# get list of Ne files output 
with open("%s/my_Ne_est_list.txt"%ddir, "w") as  f:
    for rfile in rfiles:
	for ind in range(subset_interval,num_indv,subset_interval): #the same subset as estimating Ne in previous step
	    emfile = "%s/%s.ind%d.EMprobs.out" % (ddir, rfile, ind)
	    f.write("%s %s\n" % (rfile, emfile))


with open(script, "w") as f:
    f.write("/home-fn/users/nscc1082/software/fs-2.0.7/scripts/neaverage.pl -o ./%s/my_Ne_est.txt -l ./%s/my_Ne_est_list.txt\n" % (ddir, ddir))
os.system("chmod +x %s" % script)

