#!/usr/bin/env python2

def comb_seq(dna, pattern):
	len_pattern = len(pattern)
	len_dna =  len(dna)

	index_pos = []
	for i in xrange(0,len_dna):
		if dna[i:i+len_pattern] == pattern:
			index_pos.append(str(i + 1))
	return ' '.join(index_pos)



if __name__ == '__main__':
    dataset = open("/Users/krishna/Downloads/rosalind_subs.txt").read().split("\n")    
    #dataset = 'GATATATGCATATACTT\nATAT'.split("\n")

    dna = comb_seq(dataset[0],dataset[1])
    print dna
			