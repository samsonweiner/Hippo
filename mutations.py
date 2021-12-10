import numpy as np
from tree import Node, Tree

class SNV():
    def __init__(self, name=None, chrom=None, allele=None, pos=None, ref_nuc=None, alt_nuc=None):
        self.name = name
        self.chrom = chrom
        self.allele = allele
        self.pos = pos
        self.ref_nuc = ref_nuc
        self.alt_nuc = alt_nuc

class CNV():
    def __init__(self, name=None, chrom=None, allele=None, start=None, length=None, amp=None):
        self.name = name
        self.chrom = chrom
        self.allele = allele
        self.start = start
        self.length = length
        #if amp == 0, that means deletion event
        self.amp = amp


#sets the branch length to the root as scale*total branch length in the tree.
def set_root_branchlen(tree, scale):
    tot = tree.get_total_branchlen()
    tree.root.set_len(tot*scale)

def place_SNVs_on_tree(tree, num_snv):
    tot = tree.get_total_branchlen()
    nodes = [node for node in tree.iter_descendants()]
    placements = np.random.choice(nodes, num_snv, p=[node.length/tot for node in nodes])
    return placements

def place_CNVs_on_tree(tree, num_cnv):
    tot = tree.get_total_branchlen()
    nodes = [node for node in tree.iter_descendants()]
    placements = np.random.choice(nodes, num_cnv, p=[node.length/tot for node in nodes])
    return placements

def get_chrom_proportions(chrom_len):
    total_len = sum(chrom_len.values())
    return dict(zip(chrom_len.keys(), [val/total_len for val in chrom_len.values()]))

def mutate_nucleotide(ref_nuc):
    nucleotides = ['A', 'C', 'G', 'T']
    nucleotides.remove(ref_nuc)
    alt_nuc = np.random.choice(nucleotides, 1)[0]
    return alt_nuc

def create_SNV_events(tree, num_snv, placements, snv_names, chrom_proportions, chrom_len, chrom_seq):
    for i in range(num_snv):
        node = placements[i]

        chrom = np.random.choice(list(chrom_len.keys()), 1, p=list(chrom_proportions.values()))[0]
        allele = np.random.binomial(1, 0.5)
        pos = np.random.randint(chrom_len[chrom])
        ref_nuc = chrom_seq[chrom][pos]
        alt_nuc = mutate_nucleotide(ref_nuc)

        new_SNV = SNV(name=snv_names[i], chrom=chrom, allele=allele, pos=pos, ref_nuc=ref_nuc, alt_nuc=alt_nuc)
        node.SNVs.append(new_SNV)

#TODO: update chromosome lengths and sequences when adding copy numbers. Might need to do this by traversing the tree...

def create_CNV_events(tree, num_cnv, placements, cnv_names, amp_rate, del_rate, min_size, alpha, beta, chrom_proportions, chrom_len, chrom_seq):
    for i in range(num_cnv):
        node = placements[i]

        chrom = np.random.choice(list(chrom_len.keys()), 1, p=list(chrom_proportions.values()))[0]
        allele = np.random.binomial(1, 0.5)
        size = int(np.floor(np.random.gamma(alpha, beta)*1000000 + min_size))
        start = np.random.randint(chrom_len[chrom] - size)
        #CN deletion
        if np.random.binomial(1, del_rate) == 1:
            extra_copies = 0
        #CN gain
        else:
            extra_copies = np.random.geometric(amp_rate)
        
        new_CNV = CNV(name=cnv_names[i], chrom=chrom, allele=allele, start=start, length=size, amp=extra_copies)
        node.CNVs.append(new_CNV)
