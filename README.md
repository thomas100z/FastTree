# FastTree

4 key points of FastTree

## 1. storing profiles instead of distance matrices
Neighbour joining is done by storing profiles for internal nodes. 

## 2. Neighbour joining performed by two heuristcs
Decrease amount of joins considered

## 3. initial toplogy refined by Neigherst Neighbour Interchange (NNI)

## 4. Bootstrap 
???

## FastTree algorithm
Workings of FastTree in order:
1. read alignment file and construct nodes.
2. calculate total profile
3. calculate top hits for each node
4. construct initial tree topology 
5. refine initial topology by NNI
6. iterativly improve topology
7. final branch lengths
8. support values
