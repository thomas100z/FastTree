import argparse

parser = argparse.ArgumentParser(description='Re-creation of FastTree Algorithm in Python, original authors ...')
parser.add_argument('input file', metavar='input_file', type=str,
                    help='document containing nucleotide sequence in .aln format')

args = parser.parse_args()
