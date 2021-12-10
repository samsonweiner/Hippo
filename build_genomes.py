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
        self.amp = amp  #if amp == -1, that means deletion event


#sets the branch length to the root as scale*total branch length in the tree.
def set_root_branchlen(tree, scale):
    tot = tree.get_total_branchlen()
    tree.root.set_len(tot*scale)

def place_SNVs_on_tree(tree, num_snv):
    tot = tree.get_total_branchlen()
    nodes = [node for node in tree.iter_descendants()]
    placements = np.random.choice(nodes, num_snv, p=[node.length/tot for node in nodes])
    snv_map = {}
    for i, node in enumerate(placements):
        if node in snv_map:
            snv_map[node].append(i)
        else:
            snv_map[node] = [i]
    return snv_map

def place_CNVs_on_tree(tree, num_cnv):
    tot = tree.get_total_branchlen()
    nodes = [node for node in tree.iter_descendants()]
    placements = np.random.choice(nodes, num_cnv, p=[node.length/tot for node in nodes])
    cnv_map = {}
    for i, node in enumerate(placements):
        if node in cnv_map:
            cnv_map[node].append(i)
        else:
            cnv_map[node] = [i]
    return cnv_map

def get_chrom_proportions(sequence):
    combined_chrom_len = dict(zip(sequence.keys(), map(lambda x: len(sequence[x][0]) + len(sequence[x][1]), sequence.keys())))
    total_len = sum(combined_chrom_len.values())
    combined_chrom_proportions = dict(zip(combined_chrom_len.keys(), [val/total_len for val in combined_chrom_len.values()]))
    return combined_chrom_proportions

def mutate_nucleotide(ref_nuc):
    nucleotides = ['A', 'C', 'G', 'T']
    nucleotides.remove(ref_nuc)
    alt_nuc = np.random.choice(nucleotides, 1)[0]
    return alt_nuc

# Simulates an SNV event occuring at a specific node
def create_SNV_event(node, snv_name, chrom_proportions, sequence):
        chrom = np.random.choice(list(chrom_proportions.keys()), 1, p=list(chrom_proportions.values()))[0]
        allele = np.random.binomial(1, 0.5)
        pos = np.random.randint(len(sequence[chrom][allele]))
        ref_nuc = sequence[chrom][allele][pos]
        alt_nuc = mutate_nucleotide(ref_nuc)
        node.sequence[chrom][allele] = node.sequence[chrom][allele][:pos] + alt_nuc + node.sequence[chrom][allele][pos+1:]

        new_SNV = SNV(name=snv_name, chrom=chrom, allele=allele, pos=pos, ref_nuc=ref_nuc, alt_nuc=alt_nuc)
        node.SNVs.append(new_SNV)

# Simulates a CNV event occuring at a specific node
def create_CNV_event(node, cnv_name, amp_rate, del_rate, min_size, size_rate, chrom_proportions, sequence):
        chrom = np.random.choice(list(chrom_proportions.keys()), 1, p=list(chrom_proportions.values()))[0]
        allele = np.random.binomial(1, 0.5)
        size = int(np.floor(np.random.exponential(size_rate))) + min_size
        start = np.random.randint(len(sequence[chrom][allele]) - size)
        #CN deletion
        if np.random.binomial(1, del_rate) == 1:
            extra_copies = -1
        #CN gain
        else:
            extra_copies = np.random.geometric(amp_rate)
        
        if extra_copies == -1:
            node.sequence[chrom][allele] = node.sequence[chrom][allele][:start] + node.sequence[chrom][allele][start+size:]
        else:
            node.sequence[chrom][allele] = node.sequence[chrom][allele][:start] + (1+extra_copies)*node.sequence[chrom][allele][start:start+size] + node.sequence[chrom][allele][start+size:]
        
        new_CNV = CNV(name=cnv_name, chrom=chrom, allele=allele, start=start, length=size, amp=extra_copies)
        node.CNVs.append(new_CNV)


def build_genomes(tree, num_snv, num_cnv, root_scale, chrom_seq, amp_rate, del_rate, min_size, size_rate):
    set_root_branchlen(tree, root_scale)
    snv_map = place_SNVs_on_tree(tree, num_snv)
    cnv_map = place_CNVs_on_tree(tree, num_cnv)

    tree.root.init_sequence(chrom_seq)
    chrom_proportions = get_chrom_proportions(tree.root.sequence)

    for node in tree.iter_descendants():
        if not node.is_root():
            node.inheret()
        #Create cnv/snv events along the edge above the node. Update the genome sequence.
        if node in snv_map:
            chrom_proportions = get_chrom_proportions(node.sequence)
            for i in snv_map[node]:
                create_SNV_event(node, 'snv' + str(i), chrom_proportions, node.sequence)
        if node in cnv_map:
            for i in cnv_map[node]:
                chrom_proportions = get_chrom_proportions(node.sequence)
                create_CNV_event(node,'cnv' + str(i), amp_rate, del_rate, min_size, size_rate, chrom_proportions, node.sequence)
