#!/usr/bin/env python2
# Completing a Tree

if __name__ == '__main__':
	dataset = open("/Users/krishna/Downloads/rosalind_tree_y.txt").read().split('\n')
	n = int(dataset[0])
	adj_list = dataset[1:]
	edges = n - len(adj_list)
	print edges