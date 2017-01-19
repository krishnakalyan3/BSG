#!/usr/bin/env python2
# Creating a Character Table

from numpy import zeros
from Newik import *

def char_table(input_nk):
	tree_nk = Newick(input_nk)
	node_n = lambda n: 'Node' not in n
	edge_noname = lambda i: 'Node_' in i[0] and 'Node_' in i[1]
	nodeo = {name:index for index,name in enumerate(sorted(filter(node_n, [node.name for node in tree_nk.nodes])))}
	edge = filter(edge_noname, tree_nk.edge_names())
	M = zeros((len(edge), len(nodeo)), dtype=int)
	for i, edge in enumerate(edge):
    		taxa = filter(node_n, set(tree_nk.get_descendants(edge[0])) & set(tree_nk.get_descendants(edge[1])))
    		for t in taxa:
        		M[i][nodeo[t]] = 1
	return M


if __name__ == '__main__':
	#newick_in = '(dog,((elephant,mouse),robot),cat);'
	newick_in = open("/Users/krishna/Downloads/rosalind_ctbl.txt").read().strip()
	table =  char_table(newick_in)
	for i in table:
		print "".join(str(x) for x in i)

