#!/usr/bin/env python2
import itertools


def word(data, n):
    k =[]
    for i in range(1,n + 1):
         data_sub = itertools.product(data, repeat=i)
         for j in data_sub:
            k.append(''.join(j))
    
    srt_perm = sorted(k,key=lambda word: [data.index(c) for c in word])
    return srt_perm

if __name__ == '__main__':
    #dataset = "D N A\n3".split("\n")
    dataset = open("/Users/krishna/Downloads/rosalind_lexv.txt").read()
    data_split = dataset.split("\n")
    list_word = data_split[0].split(' ')
    for i in word(list_word, int(data_split[1])):
       print i

