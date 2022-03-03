import argparse
import os
import subprocess

from cell_lineage_tree import *
from process_ref import process_ref
from build_genomes import build_genomes
from sim_reads import *
from tree import Node, Tree

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--module', type=str, default='all', help='Select one of three modules.')
parser.add_argument('-r', '--ref-path', type=str, default='hg38.fa', help='Reference file in fasta format.')
parser.add_argument('-o', '--out-path', type=str, default='/', help='Path to output directory.')
parser.add_argument('-t', '--cell-tree', type=str, default='ms')
parser.add_argument('-n', '--num-cells', type=int, default=10)
parser.add_argument('-S', '--num-snv', type=int, default=50)
parser.add_argument('-C', '--num-cnv', type=int, default=10)
parser.add_argument('-e', '--min-cn-length', type=int, default=500000)
parser.add_argument('-s', '--cn-length-rate', type=float, default=0.000001)
parser.add_argument('-a', '--cn-amp-rate', type=float, default=0.5)
parser.add_argument('-d', '--cn-del-rate', type=float, default=0.5)
parser.add_argument('-b', '--root-branch', type=float, default=0.5)
parser.add_argument('-u', '--mu', type=int, default=60)
parser.add_argument('-y', '--disp', type=int, default=120)
parser.add_argument('-x', '--region-length', type=int, default=200000)
parser.add_argument('-l', '--read-length', type=int, default=35)
parser.add_argument('-w', '--wgsim-path', type=str, default='../wgsim/wgsim')
parser.add_argument('-f', '--ms-path', type=str, default='../msdir/ms')
parser.add_argument('-p', '--num-processors', type=int, default=1)
args = parser.parse_args()

module = args.module
ref_path = args.ref_path
out_path = args.out_path
tree_type = args.cell_tree
num_cells = args.num_cells
num_snv = args.num_snv
num_cnv = args.num_cnv
min_size = args.min_cn_length
size_rate = args.cn_length_rate
amp_rate = args.cn_amp_rate
del_rate = args.cn_del_rate
root_scale = args.root_branch
mu = args.mu
disp = args.disp
region_len = args.region_length
read_len = args.read_length
wgsim_path = args.wgsim_path
ms_path = args.ms_path

if not os.path.isdir(out_path):
    call = subprocess.call(['mkdir', out_path])

if module == 'tree':
    pass

if module == 'genomes':
    pass

if module == 'reads':
    pass

if module == 'all':

    #Module 1
    if tree_type == 'ms':
        tree_str = call_ms(ms_path, num_cells, out_path)
        t = Tree(newick=tree_str)
        for leaf in t.iter_leaves():
            leaf.name = 'leaf' + leaf.name

    elif tree_type == 'random':
        t = gen_random_topology(num_cells)
        add_branchlen_ultrametric(t)
        add_branchlen_deviation(t)

    #Module 2
    chrom_seq = process_ref(ref_path)
    build_genomes(t, num_snv, num_cnv, root_scale, chrom_seq, amp_rate, del_rate, min_size, size_rate)

    #Module 3
    for leaf in t.iter_leaves():
        write_seq(leaf, out_path)
        seq_lens = get_seq_lens(leaf)
        simulate_reads(leaf, mu, disp, region_len, read_len, seq_lens, out_path, wgsim_path)


