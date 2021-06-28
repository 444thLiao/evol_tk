"""
Provide download function with advanced modifications.

Currently, it provides API of FTP. 
"""

from ftplib import FTP
from tqdm import tqdm
from os.path import join
links = ["ftp.sra.ebi.ac.uk/vol1/fastq/ERR341/003/ERR3412973/ERR3412973_1.fastq.gz",
"ftp.sra.ebi.ac.uk/vol1/fastq/ERR341/003/ERR3412973/ERR3412973_2.fastq.gz"]

def get_ftp_file(link,odir='./'):
    server,_,path = link.partition('/')
    path,_,filename = path.rpartition('/')
    ftp=FTP()
    ftp.connect(server) 
    ftp.login()
    ftp.cwd(path)
    
    ## below codes provide transvering and searching for targeted suffix
    # all_f = ftp.nlst()
    # protein_f = [_ for _ in all_f if _.endswith('_protein.faa.gz')]
    # if protein_f:
    #     protein_f = protein_f[0]
    #     with open(f'./{protein_f}','wb') as f1:
    #         ftp.retrbinary(f'RETR {protein_f}',f1.write,1024)
    # else: 
    #     gbk = [_ for _ in all_f if _.endswith('genomic.gbff.gz')][0]
    #     with open(f'./{gbk}','wb') as f1:
    #         ftp.retrbinary(f'RETR {gbk}',f1.write,1024)
    
    with open(join(odir,filename), 'wb') as fd:
        total = ftp.size(filename)
        with tqdm(total=total) as pbar:
            def callback_(data):
                l = len(data)
                pbar.update(l)
                fd.write(data)

            ftp.retrbinary('RETR {}'.format(filename), callback_,1024)
    # ftp.retrbinary("RETR " + filename, open(odir+ filename, 'wb').write)
    ftp.quit()            
            