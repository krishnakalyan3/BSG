#!/usr/bin/env python2

import fasta as f
import problem1 as p1
import problem13 as p13
import problem2 as p2

def remove_sub2(dna_string,patterns):
	for pat in patterns:
		dna_string = dna_string.replace(pat,'')
	return dna_string

if __name__ == '__main__':
    dataset = open("/Users/krishna/Downloads/rosalind_splc.txt").read().split('>')
    #dataset = open("/Users/krishna/MIRI/BSG/data/splice.txt").read().split('>')
    parse_data = f.parse_fasta(dataset)
    filter_data = map(lambda (x,y) : y, parse_data)
    trans_str = remove_sub2(filter_data[0], filter_data[1:])
    trans_rna_str = p1.dna_rna(trans_str)
    print p2.code_protein(trans_rna_str)