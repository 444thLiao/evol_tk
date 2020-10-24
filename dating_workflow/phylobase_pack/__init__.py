



# mcmctree usage of calibrations
# GCA_000011385.1|GCA_002746535.1 >33<35  root
# GCA_000011385.1|GCA_000013205.1 >23.2<30        gleobacteria
# GCA_000196515.1|GCA_001548455.1 >16<19  Nostocales
# GCA_000317575.1|GCA_000317025.1 >17<19  pleurocapsales

# phylobayes 
# 4
# GCA_000011385.1 GCA_003576915.1 35 33
# GCA_000011385.1 GCA_000013205.1 30 23.2
# GCA_000196515.1 GCA_001548455.1 19 16
# GCA_000317575.1 GCA_000317025.1 19 17


def convert_cal(infile):
    ori_text = open(infile).read().split('\n')
    valid_rows = [_ for _ in ori_text if _ and not _.startswith('#')]
    
    valid_rows = [_ for _ in valid_rows if '>' in _ or '<' in _]
    num_r = len(valid_rows)
    rows = [str(num_r)]
    for row in valid_rows:
        LCA = row.split('\t')[0]
        time = row.split('\t')[1].strip('>').split('<')
        n1,n2 = LCA.split('|')
        times = ['-1' if _ else _ for _ in time]
        times = times[::-1]
        rows.append(' '.join([n1,n2,]+ times))
    return '\n'.join(rows)