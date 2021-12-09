import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, required=True, help='Path to input file.')
parser.add_argument('-o', '--output', type=str, default='/', help='Path to output directory.')
parser.add_argument('-m', '--num-cells', type=int, default=20, help='Number of cells.')
args = parser.parse_args()
