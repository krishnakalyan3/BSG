#!/usr/bin/env python2
# Fixing an Inconsistent Character Set

from CBP import *

if __name__ == '__main__':
	with open('/Users/krishna/Downloads/rosalind_cset.txt') as f:
        	chars = f.read()
        	chars = chars.strip()
        	chars = chars.split()

	#chars = '100001\n000110\n111000\n100111'
	taxa = [str(n) for n in range(len(chars[0]))]
    	for i in range(len(chars)):
        	nchars = chars[:i] + chars[i+1:]
        	if consistent(nchars,taxa):
            		print("\n".join(nchars))
            		break