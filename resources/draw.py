from Bio import Phylo

tree = Phylo.read('resources/test-small.out', 'newick')
Phylo.draw(tree)