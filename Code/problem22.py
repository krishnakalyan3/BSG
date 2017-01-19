#!/usr/bin/env python2
# Finding a Shared Spliced Motif
import fasta as fs
import numpy as np

def lc_seq(d1, d2):
    matrix = np.zeros((len(d1) + 1, len(d2) + 1), dtype=int)
    for i in range(len(d1)):
        for j in range(len(d2)):
            if d1[i] == d2[j]:
                matrix[i + 1, j + 1] = matrix[i, j] + 1
            else:
                matrix[i + 1, j + 1] = max(matrix[i + 1, j], matrix[i, j + 1])
    ls = ''
    i, j = len(d1), len(d2)
    while i * j != 0:
        if matrix[i, j] == matrix[i - 1, j]:
            i -= 1
        elif matrix[i, j] == matrix[i, j - 1]:
            j -= 1
        else:
            ls = d1[i - 1] + ls
            i -= 1
            j -= 1
    return ls


if __name__ == '__main__':
	dataset = open("//Users/krishna/Downloads/rosalind_lcsq.txt").read().split('>')
	dataset_parse = fs.parse_fasta(dataset)
	print lc_seq(dataset_parse[0][1], dataset_parse[1][1])
	#print lc_seq('AACCTTGG', 'ACACTGTGA')
