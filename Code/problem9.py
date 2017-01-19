#!/usr/bin/env python2

def count_dna(data):
    nd = {}
    for i in data:
        if i not in nd:
            nd[i] = 1
        else:
            nd[i] = nd[i] + 1

    count_op = (str(nd['A']) +' '+ str(nd['C']) + ' ' + str(nd['G']) + ' ' + str(nd['T']) )
    return count_op

if __name__ == '__main__':
    dataset = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
    assert count_dna(dataset) == "20 12 17 21"

    dataset_large = open("/Users/krishna/Downloads/rosalind_dna.txt").read()
    print count_dna(dataset_large)
    
