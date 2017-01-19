#!/usr/bin/env python2

def getindices(s):
    return [i for i, c in enumerate(s) if c.isupper()]

def parse_fasta(dna):
	dna_filtred  = filter(lambda x : len(x) > 0, dna)
	
	fasts_seq = []
	for i in dna_filtred:
		dna_fix = i.replace('\n','')
		#print getindices(dna_fix)[1]
		name_label = dna_fix[0:getindices(dna_fix)[1]]
		name_dna = dna_fix[getindices(dna_fix)[1]:]
		fasts_seq.append((name_label,name_dna))
		
	return fasts_seq

if __name__ == '__main__':
    dataset = open("/Users/krishna/Downloads/rosalind_splc.txt").read().split('>')

    print parse_fasta(dataset)