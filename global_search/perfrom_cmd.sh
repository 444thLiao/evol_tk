
# amoC
# trimal
python3 ~/script/evolution_relative/global_search/build_tree_exe.py ./nr_retrieve_amoC/filtered_by_kegg.faa .iqtree.treefile
python3 ~/script/evolution_relative/global_search/get_info.py -i nr_retrieve_amoC/filtered_by_kegg.faa_aln.dir/iqtree.treefile/used_ids.list -o nr_retrieve_amoC/filtered_by_kegg.faa_aln.dir/iqtree.treefile

# amoB
# notrim 
python3 ~/script/evolution_relative/global_search/build_tree_exe.py ./nr_retrieve_amoB/filtered_by_kegg.faa .iqtree.treefile
python3 ~/script/evolution_relative/global_search/get_info.py -i nr_retrieve_amoB/filtered_by_kegg.faa_aln.dir/iqtree.treefile/used_ids.list -o nr_retrieve_amoB/filtered_by_kegg.faa_aln.dir/iqtree.treefile


# hao
python3 ~/script/evolution_relative/global_search/build_tree_exe.py ./nr_retrieve_hao/filtered_by_kegg.faa .iqtree.no_trim.treefile
python3 /home-user/thliao/script/evolution_relative/global_search/get_info.py -i ./nr_retrieve_hao/filtered_by_kegg.faa_aln.dir/iqtree.no_trim.treefile/used_ids.list -o nr_retrieve_hao/filtered_by_kegg.faa_aln.dir/iqtree.no_trim.treefile/


# amoA
python3 ~/script/evolution_relative/global_search/build_tree_exe.py ./nr_retrieve_removeENV_amoA/cluster_98 .newick
python3 /home-user/thliao/script/evolution_relative/global_search/get_info.py -i ./nr_retrieve_removeENV_amoA/cluster_98_aln.dir/newick/used_ids.list -o nr_retrieve_removeENV_amoA/cluster_98_aln.dir/newick/


# nxrA
python3 ~/script/evolution_relative/global_search/build_tree_exe.py ./nr_retrieve_nxrA/cluster_95_filtered_lengths.fa .iqtree.treefile
python3 /home-user/thliao/script/evolution_relative/global_search/get_info.py -i ./nr_retrieve_nxrA/cluster_95_filtered_lengths.fa_aln.dir/iqtree.treefile/used_ids.list -o nr_retrieve_nxrA/cluster_95_filtered_lengths.fa_aln.dir/iqtree.treefile/



python3 ~/script/evolution_relative/global_search/reannotate_tmp.py nr_retrieve_nxrA/cluster_95_filtered_lengths.fa_aln.dir/iqtree.treefile nr_retrieve_hao/filtered_by_kegg.faa_aln.dir/iqtree.no_trim.treefile/ nr_retrieve_amoB/filtered_by_kegg.faa_aln.dir/iqtree.treefile nr_retrieve_amoC/filtered_by_kegg.faa_aln.dir/iqtree.treefile