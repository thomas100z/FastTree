#!/usr/bin/env python3
from classes import *
import math
import argparse
import logging
from Bio import Phylo

# argument parser for cli
parser = argparse.ArgumentParser(description='Re-creation of FastTree Algorithm in Python, original authors ...')
parser.add_argument('input_file', metavar='input_file', type=str,
                    help='document containing nucleotide sequence in .aln format')
parser.add_argument('output_file', metavar='output_file', type=str,
                    help='output file in newick format')
parser.add_argument('-v', help='print verbose debugging information', action='store_true')
parser.add_argument('-t', help='view result in tree format', action='store_true')
args = parser.parse_args()

# logging settings
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='[%I:%M:%S %p]')
logger = logging.getLogger('FastTree')

if args.v:
    logging.basicConfig(level=logging.DEBUG)
    logger.setLevel(logging.DEBUG)

# parse the alignment file
parser = AlignmentParser(args.input_file)

# initialize Nodes and Tree
nodes = parser.get_data()
N = len(nodes)
m = round(math.sqrt(N))
tree = Tree(nodes, m, N)

a = 4
nni_round = math.log(N) / math.log(2) + 1

for node in nodes:
    logger.debug(f'{node.name=}\t{node.alignment}\t {node.profile=}')

# calculate total profile
tp = TotalProfile(nodes)

logging.debug('total profile:')
logging.debug(f'{tp.total_profile=}')

# set top hist list for every node
tree.set_top_hits()

for node in nodes:
    logging.debug(f'{node.name=} top hits: {[_.name for _ in node.top_hits]}\t {round(node.best_known.distance)=} '
                  f'{node.best_known.node.name=}')

# set the best join
tree.construct_initial_topology()

# NNI
previous_changes = 0
for i in range(round(nni_round)):
    previous_changes = tree.joins

    if tree.joins > 200:
        tp.recompute(tree.active_nodes)

    tree.nearest_neighbor_interchange()

    if previous_changes == tree.joins:
        break

# local bootstrap

# branch length
# tree.calculate_branch_length()

# print the tree
print(tree.to_newick())
tree.save(args.output_file)

if args.t:
    tree = Phylo.read(args.output_file, 'newick')
    Phylo.draw(tree)
