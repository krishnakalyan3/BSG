#!/usr/bin/env python2
# Finding a Spliced Motif
import fasta as fs


def sub_sequence(dna, sub_seq):
	index1 = 0
	index2 = 0
	res = []
	for d in dna:
		index1 += 1
		if sub_seq[index2] == d:
			res.append(str(index1))
			index2 += 1
			if index2 == len(sub_seq):
				break
	return res

if __name__ == '__main__':
	dataset = open("//Users/krishna/Downloads/rosalind_sseq.txt").read().split('>')
	dataset_parse = fs.parse_fasta(dataset)
	print ' '.join(sub_sequence(dataset_parse[0][1], dataset_parse[1][1]))
	#print sub_sequence('TATGCTAAGATC', 'ACG')