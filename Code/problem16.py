#!/usr/bin/env python2
# Finding a Protein Motif
import urllib
import re

def lcs(strings):
    strings = sorted(strings)
    short_string = strings[0]
    other_strings = strings[1:]

    len_str = len(short_string)
    m = ''
    for i in range(0, len_str):
        for j in range(len_str, i + len(m), -1):
            s1 = short_string[i:j]

            matched_all = True
            for s2 in other_strings:
                if s1 not in s2:
                    matched_all = False
                    break

            if matched_all:
                m = s1
                break 
    return m

if __name__ == '__main__':
	dataset = open("/Users/krishna/Downloads/rosalind_lcsm.txt").read().split('>')
	
	fixed_data  = []
	for i in dataset:
		data_join = ''.join(i[14:].split('\n')) 
		fixed_data.append(data_join)

	data_remove_empty = filter(lambda x : len(x)> 0, fixed_data)
	print lcs(data_remove_empty)