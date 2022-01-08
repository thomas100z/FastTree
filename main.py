#!/usr/bin/env python3
from classes import *
import math
import argparse

parser = argparse.ArgumentParser(description='Re-creation of FastTree Algorithm in Python, original authors ...')
parser.add_argument('input_file', metavar='input_file', type=str,
                    help='document containing nucleotide sequence in .aln format')
parser.add_argument('output_file', metavar='output_file', type=str,
                    help='output file in newick format')

args = parser.parse_args()

# parse the alignment file
parser = AlignmentParser(args.input_file)

# initialize Nodes and Tree
nodes = parser.get_data()
N = len(nodes)
m = round(math.sqrt(N))
tree = Tree(nodes, m, N)

a = 4
nni_round = math.log(N)/math.log(2) + 1

# calculate total profile
tp = TotalProfile(nodes)

# set top hist list for every node
tree.set_top_hits()

# set the best join
tree.construct_initial_topology()

# NNI
for i in range(round(nni_round)):

    if tree.joins > 200:
        tp.recompute()

    tree.nearest_neighbor_interchange()

# local bootstrap

# branch length
tree.calculate_branch_length()

# print the tree
print(tree.to_newick())
tree.save(args.output_file)

