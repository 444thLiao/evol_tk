from os.path import expanduser
default_name = './assembly_ids.list'
default_file = expanduser("~/.cache/ncbi-genome-download/genbank_bacteria_assembly_summary.txt")
if __name__ == '__main__':
    all_ids = open(default_name).read().split('\n')
    all_ids = set(all_ids)

    match_rows = []
    used_ids = []
    for row in open(default_file):
        _cache = row.split('\t')
        try:
            id = _cache[0]
            if id in all_ids:
                match_rows.append(row.strip('\n'))
            used_ids.append(id)
        except:
            pass
    if set(all_ids).difference(set(used_ids)):
        print('unknown id :', '\n'.join(set(all_ids).difference(set(used_ids))))
    with open('./metadata.csv','w') as f1:
        f1.write('\n'.join(match_rows))