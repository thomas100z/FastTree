#!/usr/bin/env python3
from classes import *
import math
import argparse
import logging
from Bio import Phylo

# argument parser for cli
parser = argparse.ArgumentParser(description='Re-creation of FastTree Algorithm in Python, original authors: Morgan N. '
                                             'Price, Paramvir S. Dehal, and Adam P. Arkin, FastTree: Computing Large '
                                             'Minimum Evolution Trees with Profiles instead of a Distance Matrix, '
                                             'DOI: 10.1093/molbev/msp077')

parser.add_argument('input_file', metavar='input_file', type=str,
                    help='document containing nucleotide sequence in .aln format')
parser.add_argument('output_file', metavar='output_file', type=str,
                    help='output file in newick format')
parser.add_argument('-v', help='print verbose debugging information', action='store_true')
parser.add_argument('-t', help='view result in tree format', action='store_true')

parser.add_argument('-b',metavar='bootstrap_rounds', help='Bootstrap rounds to evaluate the split of each internal node',
                    type=int, default=0)

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
L = len(nodes[0].alignment)
bootstrapping_round = args.b
tree = Tree(nodes, m, N, L, args.b, bootstrapping_round)

a = 4
nni_round = math.log(N) / math.log(2) + 1

for node in nodes:
    logger.debug(f'{node.name=}\t{node.alignment}\t {node.profile=}')

# calculate total profile
tp = TotalProfile(nodes)
tree.set_total_profile(tp)

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
for i in range(round(nni_round)):
    if tree.joins > 200:
        tp.recompute(tree.active_nodes)
        tree.set_total_profile(tp)

    tree.nearest_neighbor_interchange()

# local bootstrap
if bootstrapping_round:
    tree.bootstrap()

# branch length
tree.calculate_branch_length()

# print the tree
logger.debug(f'Final topology: {tree.to_newick()}')
tree.save(args.output_file)
print(tree.to_newick())

# show the tree if flag provided
if args.t:
    phylo_tree = Phylo.read(args.output_file, 'newick')
    Phylo.draw(phylo_tree)
