import subprocess
import os

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

    return chrom_seq

def get_chrom_lens(ref_path):
    index_path = ref_path + '.fai'
    if not os.path.isfile(index_path):
        subprocess.call(['samtools', 'faidx', ref_path])
    run = subprocess.run(['cut', '-f1,2', index_path], check=True, stdout=subprocess.PIPE, universal_newlines=True)
    lines = run.stdout.strip().split('\n')
    chrom_lens = {}
    for line in lines:
        info = line.split('\t')
        chrom_lens[info[0]] = int(info[1])
    return chrom_lens





#chrom_seq, chrom_len = process_ref('data/chr20_subset.fa')
#print(chrom_len)
#print(chrom_seq['chr20'][:200])