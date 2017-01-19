#!/usr/bin/env python2
# Counting Point Mutations

def hamming(x,y):
	count = 0
	for x1, y1 in zip(x, y):
		if (x1 != y1):
			count += 1
	return count


if __name__ == '__main__':
	dataset = open("/Users/krishna/Downloads/rosalind_hamm.txt").read()
	#dataset = 'GAGCCTACTAACGGGAT\nCATCGTAATGACGGCCT'
	split_data =  dataset.split('\n')
	print hamming(split_data[0], split_data[1])