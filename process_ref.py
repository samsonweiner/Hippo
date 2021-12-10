def process_ref(ref_path):
    chrom_seq = {}

    cur_chrom = ''
    temp_seq = ''
    with open(ref_path, 'r') as ref:
        line = ref.readline().strip('\n')
        while line:
            if '>' in line:
                if cur_chrom != '' and temp_seq != '':
                    chrom_seq[cur_chrom] = temp_seq
                cur_chrom = line[1:]
                temp_seq = ''
            else:
                temp_seq += line
            line = ref.readline().strip('\n')
    chrom_seq[cur_chrom] = temp_seq
            
    chrom_len = dict(zip(chrom_seq.keys(), map(lambda x: len(chrom_seq[x]), chrom_seq.keys())))

    return chrom_seq, chrom_len






#chrom_seq, chrom_len = process_ref('data/chr20_subset.fa')
#print(chrom_len)
#print(chrom_seq['chr20'][:200])