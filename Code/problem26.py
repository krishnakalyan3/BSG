#!/usr/bin/env python2
# Counting Unrooted Binary Trees

def count_tree(n):
	return reduce(lambda a, b: a*b % 10**6, xrange(2*n-5, 1, -2))


if __name__ == '__main__':
	dataset = open("/Users/krishna/Downloads/rosalind_cunr.txt").read()
	n = int(dataset)
	bt = count_tree(n)
	print bt
