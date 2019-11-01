
# amoC
# trimal
python3 ~/script/evolution_relative/global_search/build_tree_exe.py ./nr_retrieve_amoC/filtered_by_kegg.faa .iqtree.treefile
python3 ~/script/evolution_relative/bin/ncbi_convert/pid2all.py -i nr_retrieve_amoC/filtered_by_kegg.faa_aln.dir/iqtree.treefile/used_ids.list -o nr_retrieve_amoC/filtered_by_kegg.faa_aln.dir/iqtree.treefile/info_dir/ -f 

# amoB
# notrim 
python3 ~/script/evolution_relative/global_search/build_tree_exe.py ./nr_retrieve_amoB/filtered_by_kegg.faa .iqtree.treefile
python3 ~/script/evolution_relative/bin/ncbi_convert/pid2all.py -i nr_retrieve_amoB/filtered_by_kegg.faa_aln.dir/iqtree.treefile/used_ids.list -o nr_retrieve_amoB/filtered_by_kegg.faa_aln.dir/iqtree.treefile/info_dir/ -f 


# hao
python3 ~/script/evolution_relative/global_search/build_tree_exe.py ./nr_retrieve_hao/filtered_by_kegg.faa .iqtree.no_trim.treefile
python3 ~/script/evolution_relative/bin/ncbi_convert/pid2all.py -i ./nr_retrieve_hao/filtered_by_kegg.faa_aln.dir/iqtree.no_trim.treefile/used_ids.list -o nr_retrieve_hao/filtered_by_kegg.faa_aln.dir/iqtree.no_trim.treefile/info_dir/ -f 

# amoA(keep ENV)
python3 ~/script/evolution_relative/global_search/build_tree_exe.py ./nr_retrieve_amoA/cluster_90 .iqtree.treefile
python3 ~/script/evolution_relative/bin/ncbi_convert/pid2all.py -i ./nr_retrieve_amoA/cluster_90_aln.dir/iqtree.treefile/used_ids.list -o nr_retrieve_amoA/cluster_90_aln.dir/iqtree.treefile/info_dir/ -f 

# amoA
python3 ~/script/evolution_relative/global_search/build_tree_exe.py ./nr_retrieve_removeENV_amoA/cluster_98 .newick
python3 ~/script/evolution_relative/bin/ncbi_convert/pid2all.py -i ./nr_retrieve_removeENV_amoA/cluster_98_aln.dir/iqtree.treefile/used_ids.list -o nr_retrieve
_removeENV_amoA/cluster_98_aln.dir/iqtree.treefile/info_dir/ -f 


# nxrA
python3 ~/script/evolution_relative/global_search/build_tree_exe.py ./nr_retrieve_nxrA/cluster_95_filtered_lengths.fa .iqtree.treefile
python3 ~/script/evolution_relative/bin/ncbi_convert/pid2all.py -i ./nr_retrieve_nxrA/cluster_95_filtered_lengths.fa_aln.dir/iqtree.treefile/used_ids.list -o nr_retrieve_nxrA/cluster_95_filtered_lengths.fa_aln.dir/iqtree.treefile/info_dir/ -f 


python3 ~/script/evolution_relative/global_search/reannotate_tree.py

python3 ~/script/evolution_relative/global_search/reannotate_tmp.py nr_retrieve_nxrA/cluster_95_filtered_lengths.fa_aln.dir/iqtree.treefile nr_retrieve_hao/filtered_by_kegg.faa_aln.dir/iqtree.no_trim.treefile/ nr_retrieve_amoB/filtered_by_kegg.faa_aln.dir/iqtree.treefile nr_retrieve_amoC/filtered_by_kegg.faa_aln.dir/iqtree.treefile nr_retrieve_removeENV_amoA/cluster_98_aln.dir/iqtree.treefile