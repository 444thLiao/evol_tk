from bin.ngd import *

domain2aids, collect_info = from_name2ids("Roseobacteraceae", "refseq")
formats = ["genbank"]
odir = db_dir
downloaded_aids = []
new_domain2aids = {}
for d, aids in domain2aids.items():
    sub_aids = check_not_down(formats, aids, d, odir)
    new_domain2aids[d] = sub_aids
    downloaded_aids.extend(new_domain2aids[d])
    tqdm.write(
        f"domain: {d}, original number of ids: {len(aids)}, now ids: {len(new_domain2aids[d])} "
    )
# batch_aids = new_domain2aids['bacteria']
# ngd.download(**{"assembly_accessions": ','.join(batch_aids),
#                         "dry_run": False,
#                         "use_cache":True, # to avoid it automatic download/update the summary file
#                         "section": 'genbank',
#                         "parallel": 2,
#                         "output": db_dir,
#                         "groups": 'bacteria',  # if not assign this, it will take long time to iterate all groups
#                         "file_formats": 'genbank'})
# genbank
# original number of ids: 1795, now ids: 1149
#Thu Aug 18 16:36:39 2022


with open('./tmp.gids','w') as f1:
    f1.write('\n'.join(missing_gids))
cmd = f"python3 ~/script/evol_tk/dating_workflow/bin/postdownload.py -i /mnt/maple/thliao/data/NCBI/genbank/ -o /mnt/maple/thliao/data/NCBI/modified_data/direct_protein_files -tmp /mnt/maple/thliao/data/NCBI/modified_data/prokka_o -dry_run -gl ./tmp.gids"





