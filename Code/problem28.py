#!/usr/bin/env python2
# Counting Rooted Binary Trees

def rooted_tree(n):
	return reduce(lambda a, b: a*b % 10**6, xrange(2*n-3, 1, -2))

if __name__ == '__main__':
	dataset = open("/Users/krishna/Downloads/rosalind_root.txt").read()
	n = int(dataset)
	print rooted_tree(n)
	#print rooted_tree(4)


