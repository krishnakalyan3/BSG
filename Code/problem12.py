#!/usr/bin/env python2
from numpy import zeros

def fix_data(data):
	clean_data = []
	data_split = data.split('>')
	for i in data_split:
		clean_data.append(i.replace('\n', '')[13:])
	return filter(lambda x : len(x) >= 1, clean_data)

def matrix_gen(data):
	str_len = len(data)
	M = zeros((str_len,str_len))
	for i in xrange(0, str_len):
		for j in xrange(0, str_len):
			M[i][j] = p_distance(data[i],data[j])
	return M

def p_distance(string1,string2):
	equal = 0
	not_equal = 0
	str_len = len(string1)
	for i in xrange(0, str_len):
		if (string2[i] == string1[i]):
			equal = equal + 1
		else:
			not_equal = not_equal + 1
	return round(1.0 * not_equal/str_len,3)


if __name__ == '__main__':
	#dataset = open("/Users/krishna/MIRI/BSG/data/sample3.txt").read()
	dataset = open("/Users/krishna/Downloads/rosalind_pdst.txt").read()

	clean_data = fix_data(dataset)
	matrix = matrix_gen(clean_data)
	for i in matrix:
		print ' '.join(str(e) for e in i)