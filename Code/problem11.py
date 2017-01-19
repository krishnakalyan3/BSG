#!/usr/bin/env python2
import re
from itertools import product

def get_index(data,n):
	clean_data = []
	for i in data.split('\n')[1:]:
		clean_data.append(i)
	data_fix = ''.join(clean_data)

	kmer_list = []
	for i in xrange(0,len(data_fix)):
		kmer_list.append(data_fix[i:i+n])

	return filter(lambda x : len(x) >= n, kmer_list)

def kmer(data,n):
	kmer_dict = {}
	kmer_list =  [''.join(kmer) for kmer in list(product('ACGT', repeat = n))]
	for i in kmer_list:
		kmer_dict[i] = 0

	for i in data:
		if kmer_dict[i] == 0:
			kmer_dict[i] = 1
		else:
			kmer_dict[i] = kmer_dict[i] + 1
	return kmer_dict

if __name__ == '__main__':
  	#dataset = open("/Users/krishna/MIRI/BSG/data/sample2.txt").read()
    dataset = open("/Users/krishna/Downloads/rosalind_kmer.txt").read()
    
   	clean_data = get_index(dataset,4)
   	kmer_dict = kmer(clean_data,4)

   	comp = []
   	for key in sorted(kmer_dict):
   		comp.append(str(kmer_dict[key]))
   	print ' '.join(comp)