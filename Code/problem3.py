from os.path import dirname
import re

DNA_CODON_TABLE = {
    'TTT': 'F',     'CTT': 'L',     'ATT': 'I',     'GTT': 'V',
    'TTC': 'F',     'CTC': 'L',     'ATC': 'I',     'GTC': 'V',
    'TTA': 'L',     'CTA': 'L',     'ATA': 'I',     'GTA': 'V',
    'TTG': 'L',     'CTG': 'L',     'ATG': 'M',     'GTG': 'V',
    'TCT': 'S',     'CCT': 'P',     'ACT': 'T',     'GCT': 'A',
    'TCC': 'S',     'CCC': 'P',     'ACC': 'T',     'GCC': 'A',
    'TCA': 'S',     'CCA': 'P',     'ACA': 'T',     'GCA': 'A',
    'TCG': 'S',     'CCG': 'P',     'ACG': 'T',     'GCG': 'A',
    'TAT': 'Y',     'CAT': 'H',     'AAT': 'N',     'GAT': 'D',
    'TAC': 'Y',     'CAC': 'H',     'AAC': 'N',     'GAC': 'D',
    'TAA': '*',  	'CAA': 'Q',     'AAA': 'K',     'GAA': 'E',
    'TAG': '*',  	'CAG': 'Q',     'AAG': 'K',     'GAG': 'E',
    'TGT': 'C',     'CGT': 'R',     'AGT': 'S',     'GGT': 'G',
    'TGC': 'C',     'CGC': 'R',     'AGC': 'S',     'GGC': 'G',
    'TGA': '*',  	'CGA': 'R',     'AGA': 'R',     'GGA': 'G',
    'TGG': 'W',     'CGG': 'R',     'AGG': 'R',     'GGG': 'G'
}




def find_orf(dnadata):
	#print dnadata
	indexMstart = []
	indexstop = []

	for i in range(len(dnadata)):
		if dnadata[i] == 'M':
			indexMstart.append(i)
			

	for i in range(len(dnadata)):
		if dnadata[i] == '*':
			indexstop.append(i)
	
	# Get the last element for comparison
	#indexstop = indexstop[-1]
	# Clear list which has index more than stop / *
	#print indexMstart,indexstop
	#print indexMstart,indexstop
	#print dnadata[23]
	peptide_list = []
	seq = ""
	flag = 0
	for j in indexMstart:
		for i in range(len(dnadata)):
			if (i >= j) and (j <= indexstop):
				seq += dnadata[i]
				if (dnadata[i] == "*"):
					peptide_list.append(seq)
					seq = ""
					break	
	return peptide_list
    

	
def protein_mapping(data):
	protein = ""
	for dna in data:
		protein += DNA_CODON_TABLE[dna]
	return protein

def get_sequence(string):
	reverse = reverse_comp(string)
	forward = string

	# Forward Sequence
	codons_seq = []
	for i in [0,1,2]:
		codons_seq.append(forward[i:])

	# Reverse Sequence
	for i in [0,1,2]:
		codons_seq.append(reverse[i:])

	return codons_seq

def reverse_comp(dnadata):
	comlpement = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
	rev = dnadata[::-1]
	revcomp = ''.join(map(lambda x : comlpement[x], rev))
	return revcomp

def get_traid_mapping(dna):
	splitdna = re.findall('...',dna)
	protein = ""
	for dna in splitdna:
		protein += DNA_CODON_TABLE[dna] 
	return protein

def flattern(A):
    rt = []
    for i in A:
        if isinstance(i,list): rt.extend(flattern(i))
        else: rt.append(i)
    return rt

if __name__ == '__main__':
	small_data = "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"
	dataset = open(dirname(__file__) + "/data/rosalind_orf.txt").read()

	codon_seq = get_sequence(dataset)
	triad = []
	for codon in codon_seq:
		triad.append(get_traid_mapping(codon))

	seq_ans = []
	for i in triad:
		seq_ans.append(find_orf(i))

	dedup_seq = set(flattern(seq_ans))
	remove_stop = []
	for i in dedup_seq:
		print i[:-1]

