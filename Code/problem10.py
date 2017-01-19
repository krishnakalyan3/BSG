#!/usr/bin/env python2
import problem9 as p9

def gc_value(data):
	count_dna = map(lambda x : int(x), p9.count_dna(data).split(" "))
	op = round(float(count_dna[1]+ count_dna[2])/sum(count_dna),9) * 100
	return (op)

def data_parser(data):
	ip_sample = {}
	data_split =  data.split(">")	

	for i in filter(len, data_split):
		dna = i.replace('\n', '')
		ip_sample[dna[:13]] = gc_value(dna[13:])

	max_gc = max(ip_sample, key=ip_sample.get)

	return max_gc,ip_sample[max_gc] 



if __name__ == '__main__':
    dataset = open("/Users/krishna/MIRI/BSG/data/problem10.txt").read()
    #dataset_large = open("/Users/krishna/Downloads/rosalind_gc.txt").read()

    for i in data_parser(dataset):
    	print i

