

# cd "./dating_for/calibrations_set"

def process_time(time):
    min_v = "-1"
    max_v = "-1"
    if time.startswith('>'):
        max_v = time.split('<')[0].strip('>')
    if '<' in time:
        min_v = time.split('<')[-1].strip("><.")
    return max_v,min_v

def transform_cal_files(f):
    # from previous defined calibration file into ..... phylobayes acceptable file
    rows = [_ for _ in open(f).read().split('\n') if _]
    text = f'{len(rows)}\n'
    for row in rows:
        c_rows = row.split('\t')
        nodes = c_rows[0].split('|')
        time = c_rows[1]
        min_v, max_v = process_time(time)
        text += ' '.join(nodes) + f' {max_v} {min_v}\n'
    return text