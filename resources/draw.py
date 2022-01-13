from Bio import Phylo

tree = Phylo.read('resources/adjusted.out', 'newick')
Phylo.draw(tree)