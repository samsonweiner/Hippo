import numpy as np
import subprocess

def gen_region_coverage(mu, disp):
    prob = disp / (disp + mu)
    shape = disp
    scale = (1 - prob) / prob
    mean = np.random.gamma(shape, scale)
    depth = np.random.poisson(mean)
    return depth

def compute_coverage(depth, region_length, read_length):
    return (depth * 2 * read_length) / region_length

def get_seq_lens(leaf):
    seq_lens = {}
    for chrom in leaf.sequence:
        seq_lens[chrom] = {0: len(leaf.sequence[chrom][0]), 1: len(leaf.sequence[chrom][1])}
    return seq_lens

def write_seq(leaf, out_path):
    f0 = open(out_path + '/' + leaf.name + '_0_genome.fa', 'w+')
    f1 = open(out_path + '/' + leaf.name + '_1_genome.fa', 'w+')
    for chrom in leaf.sequence:
        f0.write('>' + chrom + '\n')
        cur_pos = 0
        tot = len(leaf.sequence[chrom][0])
        while cur_pos < tot:
            f0.write(leaf.sequence[chrom][0][cur_pos:cur_pos+50] + '\n')
            cur_pos += 50

        f1.write('>' + chrom + '\n')
        cur_pos = 0
        tot = len(leaf.sequence[chrom][1])
        while cur_pos < tot:
            f1.write(leaf.sequence[chrom][1][cur_pos:cur_pos+50] + '\n')
            cur_pos += 50
        
def simulate_reads(leaf, mu, disp, region_len, read_len, seq_lens, out_path, wgsim_path):
    for allele in [0, 1]:
        #prepare sequence for reads sim
        seq_path = out_path + '/' + leaf.name + '_' + str(allele) + '_genome.fa'
        run = subprocess.run(['samtools', 'faidx', seq_path])

        whole_forward = out_path + '/' + leaf.name + '_allele' + str(allele) + '_forward.fq'
        whole_backward = out_path + '/' + leaf.name + '_allele' + str(allele) + '_backward.fq'
        forw_file = open(whole_forward, 'a+')
        back_file = open(whole_backward, 'a+')

        for chrom in seq_lens:
            #chrom_forward = out_path + '/c-' + leaf.name + '_allele' + str(allele) + '_forward.fq'
            #chrom_backward = out_path + '/c-' + leaf.name + '_allele' + str(allele) + '_backward.fq'

            #divide into regions and compute the coverage
            total_len = seq_lens[chrom][allele]
            start = 0
            while start < total_len:
                end = start + region_len - 1
                if end > total_len:
                    end = total_len
                depth = gen_region_coverage(mu, disp)
 
                #setup sequence of region alone
                region_fa_path = out_path + '/' + leaf.name + '_allele' + str(allele) + '_region' + str(start) + '.fa'
                f = open(region_fa_path, 'w+')
                run = subprocess.run(['samtools', 'faidx', seq_path, chrom + ':' + str(start) + '-' + str(end)], stdout=f)

                #call wgsim
                temp_forward = out_path + '/temp_forward.fq'
                temp_backward = out_path + '/temp_backward.fq'
                run = subprocess.run(['./' + wgsim_path, '-h', '-N', str(depth), '-1', str(read_len), '-2', str(read_len), region_fa_path, temp_forward, temp_backward])

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

        run = subprocess.run(['rm', seq_path])
        run = subprocess.run(['rm', seq_path + '.fai'])

        
