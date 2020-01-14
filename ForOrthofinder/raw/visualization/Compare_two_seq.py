import glob,os
def visulation_loss(first_batch,storge=None):
    #todo fix all fixed variable.

    # first_batch = ['145', '142']
    if not storge:
        storge = '/tmp/temp_0626'
        if not os.path.isdir(storge):
            os.makedirs(storge)
    command = ''
    gbk74 = '/home/liaoth/data2/project/NR-Pfam/protein_db_blast/ref_gene.fasta'
    gbk = '/home/liaoth/project/NR-Pfam/New_data/extracted_part/{sample}/{sample}_*.gbk'

     #   '/home/liaoth/project/NR-Pfam/protein_db_blast/each_genome_core_gene_neighbour_gbk2/{sample}/{sample}_core_gene.gbk'
    tab = '{storge}/{sample2}_{sample1}_blastn.tab'

    # 74 is the smegmatic mc2 155 genome.
    for _i, _v in enumerate(first_batch[:-1]):
        first = _v
        second = first_batch[_i + 1]

        if not glob.glob(gbk.format(sample=first)):
            gbk_b = '/home/liaoth/project/NR-Pfam/New_data/extracted_part/{sample}/backup/{sample}_*.gbk'
            cache_gbk1s = glob.glob(gbk_b.format(sample=first))
        else:
            cache_gbk1s = [glob.glob(gbk.format(sample=first))[0]]
        cache_gbk2 = glob.glob(gbk.format(sample=second))[0]
        for74_tab = tab.format(sample2=second, sample1='74', storge=storge)
        if not os.path.isfile(for74_tab):
            dbid = glob.glob('/home/liaoth/project/NR-Pfam/New_data/extracted_part_inter_blastn/blast_db/{dbid}*.fasta*'.format(dbid=second))[0]
            cmdline = 'blastn -task blastn -query /home/liaoth/data2/project/NR-Pfam/protein_db_blast/ref_gene.fasta -db {db} -outfmt 6 -out {storge}/{outtab}'.format(db=dbid.rpartition('.')[0],storge=storge,outtab='{t1}_74_blastn.tab'.format(t1=second))
            os.system(cmdline)

        for cache_gbk1 in cache_gbk1s:
            cache_tab = tab.format(sample2=first, sample1=second,storge=storge)
            dbid = glob.glob('/home/liaoth/project/NR-Pfam/New_data/extracted_part_inter_blastn/blast_db/{dbid}*.fasta*'.format(dbid=second))[0]
            cmdline = 'blastn -query {fasta1} -db {db} -outfmt 6 -out {storge}/{outtab} -task blastn'.format(db=dbid.rpartition('.')[0],fasta1=cache_gbk1.replace('.gbk','.fasta'),storge=storge,outtab='{t1}_{t2}_blastn.tab'.format(t1=first,t2=second))
            # print(cmdline)
            os.system(cmdline)

            command = ' ' + ' '.join([cache_gbk1,cache_tab,cache_gbk2,for74_tab,gbk74])
            print('java -jar /home/liaoth/tools/act/act.jar ' + command)

            os.system('java -jar /home/liaoth/tools/act/act.jar ' + command)
            #import pdb;pdb.set_trace()

visulation_loss(['137','142'])