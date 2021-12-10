from tree import *
from mutations import *
from process_ref import *
import numpy as np

seg = 'ACGT'

seg2 = 2*seg

txt = ",,,,,rrttgg.....banana....rrr"

x = txt.strip(",.grt")

#print(x)


test_str2 = "((a,(b,d)n3)n2,c)n1;"
test_str3 = "((15:0.060,(16:0.011,20:0.011):0.049):1.035,((10:0.037,(4:0.016,18:0.016):0.020):0.751,((((11:0.022,(17:0.020,(1:0.001,13:0.001):0.019):0.001):0.052,(5:0.031,(14:0.009,19:0.009):0.023):0.042):0.091,(8:0.092,12:0.092):0.073):0.073,((6:0.021,(3:0.013,9:0.013):0.008):0.147,(2:0.030,7:0.030):0.138):0.069):0.550):0.307);"
t = Tree(newick=test_str3)

#set_root_branchlen(t, 0.5)

#place_SNVs_on_tree(t, 50)

chrom_seq, chrom_len = process_ref('data/test.fa')
chrom_proportions = get_chrom_proportions(chrom_len)

placements = place_SNVs_on_tree(t, 50)
snv_names = ['snv' + str(i) for i in range(50)]


create_SNV_events(t, 50, placements, snv_names, chrom_proportions, chrom_len, chrom_seq)

print(np.floor(2033978.8008713396)))