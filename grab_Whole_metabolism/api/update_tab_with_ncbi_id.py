from Bio import Entrez
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm

Entrez.email = 'l0404th@gmail.com'

a = pd.read_excel('manually_curated_N_cycle_genes.xlsx')
for _, row in tqdm(a.iterrows(), total=a.shape[0]):
    aa_seq = row['AA sequence']
    if not pd.isna(row['AA sequence(seq)']):
        # got before
        continue
    if pd.isna(aa_seq):
        aa_seq = row['nucl sequence']
    if len(aa_seq) <= 100:
        pid = aa_seq
        record = Entrez.read(Entrez.esearch(term=pid, db='protein'))
        if record.get('IdList', []):
            id = ','.join([str(_) for _ in record.get('IdList', [])])
            protein_seq = Entrez.efetch(db='protein', id=id, retmode='text', rettype='fasta')
            protein_seq = SeqIO.read(protein_seq, format='fasta')
            a.loc[_, 'AA sequence(seq)'] = str(protein_seq.seq)
            nuc_record = Entrez.read(Entrez.elink(dbfrom='protein', id=id, db='nuccore'))
            if len(nuc_record) != 0:
                nuc_record = nuc_record[0]
                _cache = nuc_record.get('LinkSetDb', [{}])
                if not _cache:
                    continue
                nuc_id = nuc_record.get('LinkSetDb', [{}])[0].get('Link', [{}])[0].get('Id')
                nuc_seq = Entrez.efetch(db='nuccore', id=nuc_id, retmode='text', rettype='fasta')
                nuc_seq = SeqIO.read(nuc_seq, format='fasta')
                a.loc[_, 'nucl sequence(seq)'] = str(nuc_seq.seq)
            else:
                continue
        else:
            record = Entrez.read(Entrez.esearch(term=pid, db='nuccore'))
            if record.get('IdList', []):
                id = ','.join([str(_) for _ in record.get('IdList', [])])
                nuc_seq = Entrez.efetch(db='nuccore', id=id, retmode='text', rettype='gb')
                nuc_seq = SeqIO.read(nuc_seq, format='genbank')
                a.loc[_, 'nucl sequence(seq)'] = str(nuc_seq.seq)
                pro_record = Entrez.read(Entrez.elink(dbfrom='nuccore', id=id, db='protein'))
                if len(pro_record) != 0:
                    protein_record = pro_record[0]
                    _cache = protein_record.get('LinkSetDb', [{}])
                    if not _cache:
                        continue
                    protein_id = protein_record.get('LinkSetDb', [{}])[0].get('Link', [{}])[0].get('Id')
                    protein_seq = Entrez.efetch(db='protein', id=protein_id, retmode='text', rettype='fasta')
                    protein_seq = SeqIO.read(protein_seq, format='fasta')
                    a.loc[_, 'AA sequence(seq)'] = str(protein_seq.seq)
                else:
                    continue
            else:
                continue
    else:
        a.loc[_, 'AA sequence(seq)'] = str(aa_seq)

if __name__ == '__main__':
    pass
