#!/usr/bin/env python2
# Counting Point Mutations

from numpy import zeros
import fasta as fs

def edit_distance(n1,n2):
	M = zeros((len(n1)+1,len(n2)+1), dtype=int)
	for i in range(1,len(n1)+1):
		M[i][0]= i
	for i in range(1,len(n2)+1):
		M[0][i]= i
	for i in xrange(1,len(n1)+1):
		for j in xrange(1,len(n2)+1):
			if n1[i-1] == n2[j-1]:
				M[i][j] = M[i-1][j-1]
			else:
				M[i][j] = min(M[i-1][j]+1,M[i][j-1]+1, M[i-1][j-1]+1)
	return M


if __name__ == '__main__':
	dataset = open("/Users/krishna/Downloads/rosalind_edit.txt").read().split('>')
	dataset_parse = fs.parse_fasta(dataset)
	edm = edit_distance(dataset_parse[0][1], dataset_parse[1][1])
	#print edm
	#print edm[len(dataset_parse[0][1])][len(dataset_parse[1][1])]
	print edit_distance('PLEASANTLY','MEANLY')

