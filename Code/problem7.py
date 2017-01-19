#!/usr/bin/env python2
import itertools
from os.path import dirname

def word_combination(word, n):
	word_split = word.split(" ")
	comb = itertools.product(sorted(word_split), repeat=n)
	merge = []
	for i in comb:
		merge.append(''.join(i))
	return merge
	

if __name__ == '__main__':

    dataset = open("/Users/krishna/Downloads/rosalind_lexf.txt").read()
    data = dataset.split("\n")
    comb_data = word_combination(data[0], int(data[1]))
    for i in comb_data:
    	print i