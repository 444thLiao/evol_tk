#!/home-fn/users/nscc1082/software/software/Python-2.7.9/bin/python

import os
import sys

#  perform EM inference on a subset of the data (part). We will perform 10 EM steps

min_contig_len = 10000
subset_interval = 5 # interval of subset size when perform EM inference of Ne

ddir = "parts_Ne"
if os.path.isdir(ddir):
    os.system("rm %s/*" % ddir)
else:
    os.system("mkdir %s" % ddir)


# estimate Ne
pfiles = []
num_indv = 0 # number of indivduals
for ffile in os.listdir("."):
    if ffile.endswith(".phase") and ffile.startswith("bialle_SNP.ref_"):
	with open(ffile) as f:
	    num_indv = int(f.readline()) # all files should have the same number of individuals
	    num_SNP = int(f.readline())
	    list_pos = f.readline()[:-1]
   	    # modify current phase files to chromopaiter format phase file in ./phase_for_chromopainter/
	    with open("%s/%s" % (ddir, ffile), "w") as f2:
		    list_S = "S"*num_SNP
		    f2.write("0\n%d\n%d\n%s\n%s\n%s" % (num_indv, num_SNP, list_pos, list_S, f.read()))
	    pfiles.append(ffile)

# estimate Ne for each selected contig
script = "STEP2.s2_2_part_estimate_Ne.sh"
with open(script, "w") as  f:
    for pfile in pfiles:
	rfile = pfile.replace(".phase", ".recfile")
	for ind in range(subset_interval,(num_indv+1),subset_interval):
	    f.write("echo \>\>\>IND %d of %d in %s\n" % (ind, num_indv, rfile))
            f.write("/home-fn/users/nscc1082/software/orderedPainting-master/chromopainter-0.0.4/chromopainter")
	    f.write(" -a %d %d -i 10 -in -j -g ./%s/%s -r %s -o ./%s/%s.ind%d\n\n" % (ind, ind, ddir, pfile, rfile, ddir, rfile, ind))
os.system("chmod +x %s" % script)
