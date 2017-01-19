#!/usr/bin/env python2
# Character-Based Phylogeny

from CBP import *
from nk import *

if __name__ == '__main__':
    internals = []
    
    with open('/Users/krishna/Downloads/rosalind_chbp.txt') as f:
        line = f.readline()
        line = line.strip()
        taxa = line.split()
    
        initialize(taxa)
    
        chars = []
    
        char = f.readline()
        char = char.strip()
        
        while char != "":
            split(char,taxa)
            chars.append(char)
    
            char = f.readline()
            char = char.strip()

    r = internals.pop()
    output = r.fold() + ";"

    t = newick_parse(output)
    r = edge_splits(t,taxa)
    
    for char in chars:
        assert(char in r or invert(char) in r)

    print output