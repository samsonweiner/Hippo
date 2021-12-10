from tree import *
from build_genomes import *
from process_ref import *
import numpy as np

seg = 'ACGT'

seg2 = 2*seg

txt = ",,,,,rrttgg.....banana....rrr"

x = txt.strip(",.grt")

#print(x)

test_str2 = "((a:0.687,(b:0.386,d:0.148)n3:0.095)n2:0.053,c:0.301)n1;"
test_str3 = "((15:0.060,(16:0.011,20:0.011):0.049):1.035,((10:0.037,(4:0.016,18:0.016):0.020):0.751,((((11:0.022,(17:0.020,(1:0.001,13:0.001):0.019):0.001):0.052,(5:0.031,(14:0.009,19:0.009):0.023):0.042):0.091,(8:0.092,12:0.092):0.073):0.073,((6:0.021,(3:0.013,9:0.013):0.008):0.147,(2:0.030,7:0.030):0.138):0.069):0.550):0.307);"

if False:
    t = Tree(newick=test_str2)

    #place_SNVs_on_tree(t, 50)

    #chrom_seq = process_ref('data/chr20_subset.fa')

    #build_genomes(tree, num_snv, num_cnv, root_scale, chrom_seq, amp_rate, del_rate, min_size, size_rate)
    #build_genomes(t, 100, 20, 0.5, chrom_seq, 0.5, 0.5, 10000, 50000)

    #a = [leaf for leaf in t.iter_leaves()][0]
    #print(a.SNVs)
    #print(a.CNVs)
    #print(a.sequence['chr20'][0][:1000])

disp = 120
mu = 60

prob = disp / (disp + mu)
shape = disp
scale = (1 - prob) / prob
mean = np.random.gamma(shape, scale)
depth = np.random.poisson(mean)

print(depth)