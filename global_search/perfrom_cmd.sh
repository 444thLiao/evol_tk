

python3 ~/script/evolution_relative/global_search/build_tree_exe.py ./nr_retrieve_amoC/filtered_by_kegg.faa .iqtree.treefile

python3 ~/script/evolution_relative/global_search/get_info.py -i nr_retrieve_amoC/filtered_by_kegg.faa_aln.dir/iqtree.treefile/used_ids.list -o nr_retrieve_amoB/filtered_by_kegg.faa_aln.dir/iqtree.treefile


python3 ~/script/evolution_relative/global_search/build_tree_exe.py ./nr_retrieve_hao/filtered_by_kegg.faa .iqtree.no_trim.treefile