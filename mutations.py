import numpy as np

class SNV():
    def __init__(self, name=None, chrom=None, alle=None, pos=None, ref_nuc=None, alt_nuc=None):
        self.name = name
        self.chrom = chrom
        self.alle = alle
        self.pos = pos
        self.ref_nuc = ref_nuc
        self.alt_nuc = alt_nuc

class CNV():
    def __init__(self, name=None, chrom=None, alle=None, start=None, end=None, dup_num=None):
        self.name = name
        self.chrom = chrom
        self.alle = alle
        self.start = start
        self.end = end
        self.dup_num = dup_num

