import os
from os.path import exists, join, basename
from glob import glob
from tqdm import tqdm
from collections import defaultdict
from glob import glob
import os
from os.path import *
from subprocess import check_call
from Bio import SeqIO
from tqdm import tqdm 
import pandas as pd
from ete3 import Tree
from Bio.Seq import Seq
from dating_workflow.step_script.quick_sampling import *

HOME = os.environ.get('HOME')
# database are not available for all users
mapped_mito = f"{HOME}/db/prot_db"
 #copy files from  "/home-user/sswang/project/Mito/results/euk_tree/date/acc/plant-nobc/pep-sing/"
mapped_mito_nucl = f"{HOME}/db/nucl_db"

nuc_db_dir = f"{HOME}/db/nucleus_prot_db"

dataset2usage = {"Oryza sativa": "both",
                 "Arabidopsis thaliana": "both",
                 "Physcomitrella patens": "both",
                 "Ostreococcus tauri": "both",
                 "Chondrus crispus": "both",
                 "Porphyra umbilicalis": "both",
                 "Cyanidioschyzon merolae": "both",
                 "Cyanophora paradoxa": "both",
                 "Phytophthora infestans": "both",
                 "Andalucia godoyi": "both",
                 "Jakoba libera": "mito",
                 "Jakoba bahamiensis": "mito",
                 "Reclinomonas americana": "mito",
                 "Histiona aroides": "mito",
                 "Seculamonas ecuadoriensis": "mito",
                 "Malawimonas jakobiformis": "mito",
                 "Thalassiosira pseudonana": "nuclear",
                 "Symbiodinium minutum": "nuclear",
                 "Paramecium tetraurelia": "nuclear",
                 "Oxytricha trifallax": "nuclear",
                 "Reticulomyxa filosa": "nuclear",
                 "Elphidium margaritaceum": "nuclear",
                 "Homo sapiens": "nuclear",
                 "Gallus gallus": "nuclear",
                 "Branchiostoma floridae": "nuclear",
                 "Amphimedon queenslandica": "nuclear",
                 "Candida albicans": "nuclear",
                 "Ustilago maydis": "nuclear",
                 "Pleurotus ostreatus": "nuclear",
                 "Spizellomyces punctatus": "nuclear",
                 "Dictyostelium discoideum": "nuclear",
                 "Polysphondylium pallidum": "nuclear",
                 "Acanthamoeba castellanii": "nuclear", }

color_euk = '##00E676'

### annotation mito and nuclear for bac
def anno_bac(all_ids,_outdir='./dating/mito/anno/{dataset}/anno_tab'):
    in_p_dir = "/mnt/maple/thliao/data/NCBI/modified_data/direct_protein_files"
    nuc_db_dir = nuc_db_dir
    mito_db_dir = mapped_mito
    db_dict = {'nuclear_datasets': nuc_db_dir,
               'mitoCOG': mito_db_dir}
    cmds = []
    for dataset in ['nuclear_datasets',
                    'mitoCOG', ]:
        #_outdir = f'./dating/mito/anno/{dataset}/anno_tab'
        outdir = _outdir.format(dataset=dataset)
        if not exists(outdir):
            os.makedirs(outdir)
        for aid in tqdm(all_ids):
            for db in glob(join(db_dict[dataset], '*.phr')):
                db = db.replace('.phr', '')
                in_p = join(in_p_dir, f"{aid}.faa")
                ofile = f"{outdir}/{aid}_{basename(db)}.tab"
                cmd = f"blastp -query {in_p} -db {db} -outfmt 6 -out {ofile} -evalue 1e-10 -max_hsps 1 -num_threads 20"
                # diamond_
                # cmd = f"diamond blastp -q {in_p} --quiet --db {db}.dmnd -o {ofile} -k 10 -e 1e-10 --max-hsps 1 -p 30 "
                if not exists(ofile) or getsize(ofile)==0:
                    cmds.append(cmd)

    return cmds

def parse_anno(_outdir,locus2len,name='anno_mito'):
    mi2genome2list_locus = defaultdict(dict)
    mi2genome2locus = defaultdict(dict)
    genome2mito = defaultdict(list)
    for tab in tqdm(glob(join(_outdir,'*_merged.tab'))):
        aid = basename(tab).split('_merged')[0]
        #mitoCOG = basename(tab).split('-')[-1].replace('.tab','')
        if os.path.getsize(tab)!=0:
            _df = pd.read_csv(tab,sep='\t',header=None)
            _df = _df.groupby(0).head(1)
            _df.loc[:,'gene'] = [_.split('_')[0] for _ in _df[1]]
            _df.loc[:,'align len'] = abs(_df[7] - _df[6])
            _df.loc[:,'total len'] = [locus2len.get(_,0) for _ in _df[0]]
            _df.loc[:,'cov'] = _df.loc[:,'align len']/_df.loc[:,'total len'] *100
            _df = _df.loc[(_df[10]<=1e-20) & (_df['cov']>=80),:]   # add coverage over 80
            _df = _df.sort_values(10)  # ascending
            # _df = _df.loc[_df[1].isin(mito_usage),:]  # it could restrict your blast result using genes from used euk genomes.
            gb = _df.groupby('gene')
            for gene,locus_list in gb.groups.items():
                locus_list = list(_df.loc[locus_list,0])
                all_locus = list(set(locus_list))
                mi2genome2list_locus[gene][aid] = ';'.join(all_locus)
            sub_df = gb.head(1)
            
            gene2locus = dict(zip(sub_df['gene'],sub_df[0]))
            for gene,locus in gene2locus.items():
                mi2genome2locus[gene][aid] = locus
                genome2mito[aid].append(gene)

    _df = pd.DataFrame.from_dict(mi2genome2list_locus).fillna('-')
    _bin_df = _df.applymap(lambda x: len([_ for _ in x.split(';') if _!='-']))
    with pd.ExcelWriter(join(_outdir,'..',f'{name}.xlsx')) as f1:
        _df.to_excel(f1,'locus info',index=1,index_label='Assembly IDs')
        _bin_df.to_excel(f1,'number of locus',index=1,index_label='Assembly IDs')  
    return   genome2mito,mi2genome2locus

def get_mito_prot_in_euk(mapped_mito=mapped_mito):
    mito_genes_fs = [_ for _ in glob(join(mapped_mito,'*.fas'))  ]
    used_genome2mito = defaultdict(list)
    used_mi2genome2seq = defaultdict(dict)
    for f in mito_genes_fs:
        records = list(SeqIO.parse(f,'fasta'))
        mi = basename(f).split('-')[-1].replace('.fas','')
        for _ in records:
            used_genome2mito[_.id].append(mi)
            used_mi2genome2seq[mi][_.id] = _
    return used_genome2mito,used_mi2genome2seq

def get_nuc_prot_in_euk(target_dir=nuc_db_dir):
    return get_mito_prot_in_euk(target_dir)

def remove_3rd(record,return_seq=False):
    # inplace remove the 3rd
    seq = str(record.seq)
    two_partitions = [seq[::3],seq[1::3]]
    final_seq = ''.join([''.join(_) for _ in zip(*two_partitions)])
    #final_seq = final_seq[:-2]
    if return_seq:
        return final_seq
    record.seq = Seq(final_seq)  # :-2 mean remove the lefted stop codon 
    return record

def get_mito_nucl_in_euk(mito_usage):
    mito_genes = [basename(_).split('.')[0] 
              for _ in glob(join(mapped_mito,'*.fas'))  ]
    used_mi2genome2seq_nucl = {}
    for gene in tqdm(mito_genes):
        gene_file = f'{mapped_mito_nucl}/{gene}.fas'
        # cp from /home-user/sswang/project/Mito/results/euk_tree/date/acc-cds/plant-nobc/pep-sing
        euk_genome2gene  = {_.id:_ 
                            for _ in SeqIO.parse(gene_file,'fasta') 
                            if _.id in mito_usage}
        gene = gene.split('-')[-1]
        used_mi2genome2seq_nucl[gene]= euk_genome2gene
    # find the original seq (without trimmed 3rd)        
    complete_g2genome2nucl = {}
    for gene_name,_d in used_mi2genome2seq_nucl.items():
        tmp_d = {}
        for genome,gene in _d.items():
            infaa = f'/home-user/sswang/project/Mito/data/Pisani-cds_retrieval/euk-mito/cds-query/{genome}/{gene_name}.fas'
            if not exists(infaa):
                print(f"not found {genome}  {gene_name}")
                continue
            gene_r = SeqIO.read(infaa,'fasta')
            #genome_records = list(SeqIO.parse(f'/home-user/sswang/project/Mito/data/download-mito/cds/{genome}.gene','fasta'))
            #for r in genome_records:
            if len(remove_3rd(gene_r,return_seq=1)) == len(gene.seq):
                tmp_d[genome]=gene_r
                gene_r.id = genome
            else:
                print(len(remove_3rd(gene_r,return_seq=1)),len(gene.seq))
        complete_g2genome2nucl[gene_name] = tmp_d        
    return complete_g2genome2nucl

def get_nuc_nucl_in_euk(nuclear_usage):
    euk_genome_gene2seqs = {}
    for gene_f in glob("/home-user/sswang/project/Mito/results/euk_tree/date/acc-cds/M9P-nobc/pep-sing/*.fas"):
        gene = gene_f.split('/')[-1].replace('.fas','')
        euk_genome_gene2seqs[gene] = [_ for _ in list(SeqIO.parse(gene_f,'fasta')) if _.id in nuclear_usage]
    # find the original seq (without trimmed 3rd)        
    complete_g2genome2nucl = defaultdict(list)
    for gene_name,records_list in euk_genome_gene2seqs.items():
        for gene in records_list:
            genome = gene.id
            r = SeqIO.read(f'/home-user/sswang/project/Mito/data/cds_retrieval/all/cds-query/{genome}/{gene_name}.fas','fasta')
            if remove_3rd(r,return_seq=1) == str(gene.seq):
                complete_g2genome2nucl[gene_name].append(r)
            else:
                print(gene_name,genome)      
    return complete_g2genome2nucl

