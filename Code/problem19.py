#!/usr/bin/env python2

import fasta as fs
import problem18 as ed

def algn(d1, d2, mat):
    i = len(d1)
    j = len(d2)
    t1 = []
    t2 = []
    while i > 0 and j > 0:
        if mat[i, j] == mat[i, j - 1] + 1:
            t1.insert(0, '-')
            t2.insert(0, d2[j - 1])
            j -= 1
        elif mat[i, j] == mat[i - 1, j] + 1:
            t1.insert(0, d1[i - 1])
            t2.insert(0, '-')
            i -= 1
        else:
            t1.insert(0, d1[i - 1])
            t2.insert(0, d2[j - 1])
            i -= 1
            j -= 1
    while j > 0:
        t1.insert(0, '-')
        t2.insert(0, d2[j - 1])
        j -= 1
    while i > 0:
        t1.insert(0, d1[i - 1])
        t2.insert(0, '-')
        i -= 1
    return ''.join(t1), ''.join(t2)



if __name__ == '__main__':
	dataset = open("/Users/krishna/Downloads/rosalind_edta_2.txt").read().split('>')
	dataset_parse = fs.parse_fasta(dataset)
	edm = ed.edit_distance(dataset_parse[0][1], dataset_parse[1][1])
	print edm[len(dataset_parse[0][1])][len(dataset_parse[1][1])]
	algn = algn(dataset_parse[0][1], dataset_parse[1][1], edm)
	for i in algn:
		print i
