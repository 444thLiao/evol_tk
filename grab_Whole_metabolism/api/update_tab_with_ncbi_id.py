
from Bio import Entrez
from Bio import SeqIO
Entrez.email = 'l0404th@gmail.com'

pid = 'WP_012001403'
record = Entrez.read(Entrez.esearch(term=pid, db='protein'))
if record.get('IdList', []):
    id = ','.join([str(_) for _ in record.get('IdList', [])])
    protein_seq = Entrez.efetch(db='protein',id=id,retmode='text', rettype='fasta')
    protein_seq = SeqIO.read(protein_seq,format='fasta')
    nuc_record = Entrez.read(Entrez.elink(dbfrom='protein',id=id,db='nuccore'))
    if len(nuc_record) !=0:
        nuc_record = nuc_record[0]
        nuc_id = nuc_record.get('LinkSetDb',[{}])[0].get('Link',[{}])[0].get('Id')
        nuc_seq = Entrez.efetch(db='nuccore',id=nuc_id,retmode='text', rettype='fasta')
        nuc_seq = SeqIO.read(nuc_seq,format='fasta')


t = Entrez.esearch(db='protein',term='WP_012001403')
t = Entrez.read(t)
if __name__ == '__main__':
    pass