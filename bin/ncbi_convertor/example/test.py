from bin.ncbi_convertor import NCBI_convertor
import os
os.chdir('/home-user/thliao/script/evol_tk/bin/ncbi_convertor/example')
if __name__ == "__main__":
    # test
    pids = open('./protein_ids').read().split('\n')
    convertor = NCBI_convertor(pids, db='protein')
    # convertor.check_cache(suffix=suffix, redo=redo)
    convertor.get_taxons_from_tid()
    pid2assembly_dict = convertor.pid2assembly()

    aids = open('./assembly_ids').read().split('\n')
    convertor = NCBI_convertor(aids, db='assembly')
    # convertor.check_cache(suffix=suffix, redo=redo)
    convertor.get_taxons_from_tid()

    nids = open('./nucleotide_ids').read().split('\n')
    convertor = NCBI_convertor(nids, db='nuccore')
    convertor.get_GI()
    convertor.get_db_summary()
    # convertor.check_cache(suffix=suffix, redo=redo)
    convertor.get_taxons_from_tid()
