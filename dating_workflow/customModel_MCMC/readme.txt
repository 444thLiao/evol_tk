prior to execute the following scripts.

You might want to use "/home-user/thliao/anaconda3/bin/python3.7" to avoid all potential import errors.

example: 
step1 : generating commands (pass --dryrun) OR directly submit to slurm (removing --dryrun). You could also replace the model to the other. such as LG+G.

python /home-user/thliao/script/dating_raw/customModel_MCMC.py --model 'LG+G+C60' --phy_file '/home-user/thliao/script/dating_raw/example_data/3pf.phy' --odir '~/tmp/testC60' --inbv '/home-user/thliao/script/dating_raw/example_data/in.BV' --dryrun    

step2: genearting mcmc_tree.ctl and summary in.BV from previous outputs.

python /home-user/thliao/script/dating_raw/hessian_bin.py '~/tmp/testC60' '/home-user/thliao/script/dating_raw/example_data/3pf.phy' /home-user/thliao/script/dating_raw/example_data/P39_B5P3.newick


step3: run MCMC. (optional). This script can keep running and automatically avoid the "resetting lnL" error.

python /home-user/thliao/script/dating_raw/keep_run.py ~/tmp/testC60/mcmctree/mcmc_for/