from tree import *
from build_genomes import *
from process_ref import *
from sim_reads import *
from cell_lineage_tree import *
import numpy as np
import subprocess


test_str2 = "((leafa:0.687,(leafb:0.386,leafd:0.148)n3:0.095)n2:0.053,leafc:0.301)n1;"
test_str3 = "((15:0.060,(16:0.011,20:0.011):0.049):1.035,((10:0.037,(4:0.016,18:0.016):0.020):0.751,((((11:0.022,(17:0.020,(1:0.001,13:0.001):0.019):0.001):0.052,(5:0.031,(14:0.009,19:0.009):0.023):0.042):0.091,(8:0.092,12:0.092):0.073):0.073,((6:0.021,(3:0.013,9:0.013):0.008):0.147,(2:0.030,7:0.030):0.138):0.069):0.550):0.307);"

#if False:
    #t = Tree(newick=test_str2)

    #place_SNVs_on_tree(t, 50)

    #

    #build_genomes(tree, num_snv, num_cnv, root_scale, chrom_seq, amp_rate, del_rate, min_size, size_rate)
    #build_genomes(t, 100, 20, 0.5, chrom_seq, 0.5, 0.5, 10000, 50000)

    #a = [leaf for leaf in t.iter_leaves()][0]
    #print(a.SNVs)
    #print(a.CNVs)
    #print(a.sequence['chr20'][0][:1000])

#f = open('log.txt', 'w+')
#f = open('data/test.fa', 'w+')
#run = subprocess.run(['samtools', 'faidx', 'data/chr20_subset.fa', 'chr20:0-200'], stdout=f)
#chrom_seq = process_ref('data/2chrom.fa')

#t = Tree(newick=test_str2)
#build_genomes(t, 100, 20, 0.5, chrom_seq, 0.5, 0.5, 10000, 50000)

#for leaf in t.iter_leaves():
 #   write_seq(leaf, 'data/set2')
#    seq_lens = get_seq_lens(leaf)
#    simulate_reads(leaf, 6, 12, 20000, 35, seq_lens, 'data/set2', '../wgsim/wgsim')

#call_ms('../msdir/ms', 10, 'data/t1.newick')
i = 20

f1 = open('data/chr20_full/leaf' + str(i) + '_backward.fq', 'w+')
f2 = open('data/chr20_full/leaf' + str(i) + '_forward.fq', 'w+')
call = subprocess.run(['cat', 'data/chr20_full/leaf' + str(i) + '_allele0_backward.fq', 'data/chr20_full/leaf' + str(i) + '_allele1_backward.fq'], stdout=f1)
call = subprocess.run(['cat', 'data/chr20_full/leaf' + str(i) + '_allele0_forward.fq', 'data/chr20_full/leaf' + str(i) + '_allele1_forward.fq'], stdout=f2)
call = subprocess.run(['mv', 'data/chr20_full/leaf' + str(i) + '_allele0_backward.fq', 'raw']) 
call = subprocess.run(['mv', 'data/chr20_full/leaf' + str(i) + '_allele1_backward.fq', 'raw']) 
call = subprocess.run(['mv', 'data/chr20_full/leaf' + str(i) + '_allele0_forward.fq', 'raw']) 
call = subprocess.run(['mv', 'data/chr20_full/leaf' + str(i) + '_allele1_forward.fq', 'raw']) 