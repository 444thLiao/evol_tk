import os
import sys
from os.path import dirname
sys.path.insert(0, dirname(dirname(__file__)))

# gbk_collect = []
# got_p = []
# for pfile_real in p_files_real:
#     dname = dirname(pfile_real)
#     name = basename(pfile_real).rpartition('.')[0]
#     if exists(join(dname,f'{name}.gbk')):
#         gbk_collect.append(realpath(join(dname,f'{name}.gbk')))
#         got_p.append(pfile_real)
#     elif exists(join(dirname(dname),f'{name}.gbk')):
#         gbk_collect.append(realpath(join(dirname(dname),f'{name}.gbk')))
#         got_p.append(pfile_real)
#     elif exists(join(dirname(dname),'gbkcopy',f'{name}.gbk')):
#         gbk_collect.append(realpath(join(dirname(dname),'gbkcopy',f'{name}.gbk')))
#         got_p.append(pfile_real)
#     elif exists(join(dname,f'{name}.gbk')):
#         gbk_collect.append(realpath(join(dname,f'{name}.gbk')))
#         got_p.append(pfile_real)
#     elif exists(join('/home-user/sswang/project/Rhizobiales/data/genome_source/Rhizobiales/gbk/new',f'{name}.gbk')):
#         gbk_collect.append(realpath(join('/home-user/sswang/project/Rhizobiales/data/genome_source/Rhizobiales/gbk/new',f'{name}.gbk')))
#         got_p.append(pfile_real)
#     elif exists(join('/home-user/jjtao/genome/Afipia_Pseudo_Oligo/four/',f'{name}.gbk')):
#         gbk_collect.append(realpath(join('/home-user/jjtao/genome/Afipia_Pseudo_Oligo/four/',f'{name}.gbk')))
#         got_p.append(pfile_real)
#     elif exists(join(dname,f"{name.replace('sp','sp.')}.gbk")):
#         gbk_collect.append(realpath(join(dname,f"{name.replace('sp','sp.')}.gbk")))
#         got_p.append(pfile_real)
#     else:
#         dname = dirname(dname)
#         if exists(join(dname,'gbk',f'{name}.gbk')):
#             gbk_collect.append(realpath(join(dname,'gbk',f'{name}.gbk')))
#             got_p.append(pfile_real)
# remained_ = set(p_files_real).difference(set(got_p))
# gbk_collect.append('/home-user/jjtao/genome/Bradyrhizobium/gbk/Bradyrhizobium_guangdongense_SM32.gbk')
# got_p.append('/home-user/jjtao/genome/Bradyrhizobium/gbk/Bradyrhizobium_rifense_SM32.protein')
# gbk_collect.append('/home-user/jjtao/genome/Bradyrhizobium/gbk/Bradyrhizobium_sp._R2.2-H.gbk')
# got_p.append('/home-user/jjtao/genome/Bradyrhizobium/gbk/Bradyrhizobium_sp_R2_2-H.protein')
# gbk_collect.append('/home-user/jjtao/genome/Bradyrhizobium/gbk/Bradyrhizobium_sp._cf659.gbk')
# got_p.append('/home-user/jjtao/genome/Bradyrhizobium/gbk/Bradyrhizobium_sp_cf659_CF659.protein')
# gbk_collect.append('/home-user/jjtao/genome/Rhodopseudomonas/Rhodopseudomonas_thermotolerans_JA576.gbk')
# got_p.append('/home-user/jjtao/genome/Rhodopseudomonas/Rhodopseudomonas_thermotolerans_JA576_NBRC_108863_KCTC_15144.protein')
# gbk_collect.append('/home-user/jjtao/genome/Rhodopseudomonas/Rhodopseudomonas_palustris_MAG1.gbk')
# got_p.append('/home-user/jjtao/genome/Rhodopseudomonas/Rhodopseudomonas_palustris.protein')

# for g,ori_name in zip(gbk_collect,got_p):
#     target_name = basename(ori_name).rpartition('.')[0]
#     os.system('ln -s %s /home-user/thliao/tmp_gbk_jjtao/%s.gbk' % (g,target_name))
    
# for g in glob('./gbk'):
#     cmd = f"ruby $SSW/tools/self_bao_cun/basic_process_mini/genbank2gff.rb --seq {gbk} --feature CDS > {new_file}"
    
    
# '/home-user/jjtao/backup/OrthoFinder_1/Results_Brady346/'