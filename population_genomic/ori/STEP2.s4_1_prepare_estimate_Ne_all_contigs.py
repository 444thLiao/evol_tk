#!/usr/bin/env python

import os
import sys

# based on the Ne and rm estimated via subsets of data, run chromopainter to estimate Ne again

script = "STEP2.s4_2_estimate_Ne_for_all_contigs.sh"
subset_interval = 5 # interval of subset size when perform EM inference of Ne

ddir = "parts_Ne_est"
if os.path.isdir(ddir):
    os.system("rm %s/*" % ddir)
else:
    os.system("mkdir %s" % ddir)


# use all contigs, instead of only contigs with >20000 SNPs as in step 2
pfiles = []
num_indv = 0 # number of indivduals
for ffile in os.listdir("."):
    if ffile.endswith(".phase") and ffile.startswith("bialle_SNP.ref_"):
        with open(ffile) as f:
            num_indv = int(f.readline()) # all files should have the same number of individuals
            num_SNP = int(f.readline())
            # modify current phase files to chromopaiter format phase file in ./phase_for_chromopainter/
            with open("%s/%s" % (ddir, ffile), "w") as f2:
                list_S = "S"*num_SNP
                f2.write("0\n%d\n%d\n%s\n%s\n%s" % (num_indv, num_SNP, f.readline()[:-1], list_S, f.read()))
            pfiles.append(ffile)

# load estimated Ne and rm from step 2 and step 3
globalparams = ""
with open("./parts_Ne/my_Ne_est.txt") as f:
    globalparams = f.readline()[:-1]


# estimate Ne for each contig
with open(script, "w") as  f:
    for pfile in pfiles:
        rfile = pfile.replace(".phase", ".recfile")
        for ind in range(1,(num_indv+1)):
            f.write("echo \>\>\>IND %d of %d, %s\n" % (ind, num_indv, rfile))
            f.write("/home-fn/users/nscc1082/software/orderedPainting-master/chromopainter-0.0.4/chromopainter")
            f.write(" -a %d %d %s -j -g ./%s/%s -r %s -o ./%s/%s.ind%d\n\n" % (ind, ind, globalparams, ddir, pfile, rfile, ddir, rfile, ind))
os.system("chmod +x %s" % script)
