#!/usr/bin/env python2
# Creating a Character Table from Genetic Strings

def char_str(dnal):
    tbl = set()
    for i, ch in enumerate(dnal[0]):
        array = [int(dna[i] == ch) for dna in dnal]
        if 1 < sum(array) < len(dnal)-1:
            tbl.add(''.join(map(str,array)))
    return tbl

if __name__ == '__main__':
    #dna = 'ATGCTACC\nCGTTTACC\nATTCGACC\nAGTCTCCC\nCGTCTATC'.split('\n')
    dna = open("/Users/krishna/Downloads/rosalind_cstr.txt").read().split('\n')

    table = char_str(dna)
    print '\n'.join(table)