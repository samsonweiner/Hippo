import numpy as np
import subprocess
import time

from process_ref import *

def gen_region_coverage(mu, disp):
    prob = disp / (disp + mu)
    shape = disp
    scale = (1 - prob) / prob
    mean = np.random.gamma(shape, scale)
    depth = np.random.poisson(mean)
    return depth

def write_seq(genome, name, out_path):
    f0 = open(out_path + '/' + name + '_0_genome.fa', 'w+')
    f1 = open(out_path + '/' + name + '_1_genome.fa', 'w+')
    for chrom in genome:
        f0.write('>' + chrom + '\n')
        cur_pos = 0
        tot = len(genome[chrom][0])
        while cur_pos < tot:
            f0.write(genome[chrom][0][cur_pos:cur_pos+50] + '\n')
            cur_pos += 50

        f1.write('>' + chrom + '\n')
        cur_pos = 0
        tot = len(genome[chrom][1])
        while cur_pos < tot:
            f1.write(genome[chrom][1][cur_pos:cur_pos+50] + '\n')
            cur_pos += 50
    f0.close()
    f1.close()

def simulate_reads(genome_path, name, mu, disp, region_len, read_len, seq_lens, out_path, wgsim_path):
    for allele in [0, 1]:
        #prepare sequence for reads sim
        #seq_path = out_path + '/' + name + '_' + str(allele) + '_genome.fa'
        #run = subprocess.run(['samtools', 'faidx', genome_path])

        whole_forward = out_path + '/' + name + '_allele' + str(allele) + '_forward.fq'
        whole_backward = out_path + '/' + name + '_allele' + str(allele) + '_backward.fq'
        forw_file = open(whole_forward, 'a+')
        back_file = open(whole_backward, 'a+')

        for chrom in seq_lens:

            #divide into regions and compute the coverage
            total_len = seq_lens[chrom][allele]
            #start = 125200000
            start = 0
            while start < total_len:
                stime = time.time()

                end = start + region_len - 1
                if end > total_len:
                    end = total_len
                depth = gen_region_coverage(mu, disp)

                print(depth, start, stime - time.time())
 
                #setup sequence of region alone
                region_fa_path = out_path + '/' + name + '_allele' + str(allele) + '_region' + str(start) + '.fa'
                f = open(region_fa_path, 'w+')
                run = subprocess.run(['samtools', 'faidx', genome_path, chrom + ':' + str(start) + '-' + str(end)], stdout=f)

                #call wgsim
                temp_forward = out_path + '/temp_forward.fq'
                temp_backward = out_path + '/temp_backward.fq'
                args = ' '.join([wgsim_path, '-h', '-N', str(depth), '-1', str(read_len), '-2', str(read_len), region_fa_path, temp_forward, temp_backward])
                print(args)
                run = subprocess.run([wgsim_path, '-h', '-N', str(depth), '-1', str(read_len), '-2', str(read_len), region_fa_path, temp_forward, temp_backward])


                f1 = open(temp_forward, 'r')
                forw_file.write(f1.read())
                f1.close()
                f2 = open(temp_backward, 'r')
                back_file.write(f2.read())
                f2.close()

                run = subprocess.run(['rm', temp_forward])
                run = subprocess.run(['rm', temp_backward])
                run = subprocess.run(['rm', region_fa_path])
            
                
                start = start + region_len

        #run = subprocess.run(['rm', seq_path])
        #run = subprocess.run(['rm', seq_path + '.fai'])



ref_path = '../reference/hg38.fa'
wgsim_path = 'wgsim'
ref_len = get_chrom_lens(ref_path)
for chrom in ref_len:
    chr_size = ref_len[chrom]
    ref_len[chrom] = {0: chr_size, 1: chr_size}


simulate_reads(ref_path, 'diploid1', 60, 120, 200000, 35, ref_len, '../data/WG_diploid', wgsim_path)
