import numpy as np

def gen_region_coverage(mu, disp):
    prob = disp / (disp + mu)
    shape = disp
    scale = (1 - prob) / prob
    mean = np.random.gamma(shape, scale)
    depth = np.random.poisson(mean)
    return depth

def compute_coverage(depth, region_length, read_length):
    return (depth * 2 * read_length) / region_length

