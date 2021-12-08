#!/usr/bin/env python3
import argparse
from pipeline import AlignmentParser

parser = argparse.ArgumentParser(description='Re-creation of FastTree Algorithm in Python, original authors ...')
parser.add_argument('input_file', metavar='input_file', type=str,
                    help='document containing nucleotide sequence in .aln format')

args = parser.parse_args()


sequences = AlignmentParser(args.input_file).get_data()

print("")