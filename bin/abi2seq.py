from Bio import SeqIO
import os,sys,re

def each_abi2seq(infile,ofile=None,name_pattern='(.+)__27F.*'):
    name = infile.split('/')[-1].rpartition('.')[0]
    ofile1 = f'./tmp_{name}.fastq'
    if ofile is None:
        ofile = infile.replace('.ab1','.fasta')
    SeqIO.convert(infile,'abi',ofile1,'fastq')
    
    ### modified from Biopython
    # Trims the sequence using Richard Mott's modified trimming algorithm .
    start = False  # flag for starting position of trimmed sequence
    segment = 20  # minimum sequence length
    trim_start = 0  # init start index
    cutoff = 0.001 # Q30  # default cutoff value for calculating base score
    
    seq_record = SeqIO.read(ofile1,'fastq')
    if len(seq_record) <= segment:
        seq_record = None
    else:
        # calculate base score
        score_list = [
            cutoff - (10 ** (qual / -10.0))
            for qual in seq_record.letter_annotations["phred_quality"]
        ]
        # calculate cumulative score
        # if cumulative value < 0, set it to 0
        # first value is set to 0, because of the assumption that
        # the first base will always be trimmed out
        cummul_score = [0]
        for i in range(1, len(score_list)):
            score = cummul_score[-1] + score_list[i]
            if score < 0:
                cummul_score.append(0)
            else:
                cummul_score.append(score)
                if not start:
                    # trim_start = value when cumulative score is first > 0
                    trim_start = i
                    start = True

        # trim_finish = index of highest cumulative score,
        # marking the end of sequence segment with highest cumulative score
        trim_finish = cummul_score.index(max(cummul_score))
        seq_record = seq_record[trim_start:trim_finish]
    ## Done
    
    if seq_record is not None:
        seq_record.id = re.findall(name_pattern,name)[0]
        seq_record.name = seq_record.description = ''
        SeqIO.write(seq_record,open(ofile,'w'),'fasta-2line')
    
    os.system(f"rm {ofile1} ")
    
if __name__ == '__main__':
    args = sys.argv[1:]
    
    infile = args[0]
    if len(args)==2:
        ofile = args[1]
    else:
        ofile = infile.replace('.ab1','.fasta')
    each_abi2seq(infile,ofile)
    
    from glob import glob
    if '*' in infile:
        print('There are multiple input files. ignoring output file.')
        for _infile in glob(infile):

            each_abi2seq(_infile)