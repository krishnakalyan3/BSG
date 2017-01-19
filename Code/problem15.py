#!/usr/bin/env python2
# Finding a Protein Motif
import urllib
import re

def parse_url(unipro):
	url = 'http://www.uniprot.org/uniprot/'
	fasta_ext = '.fasta'
	complete_url = url + unipro + fasta_ext
	f = urllib.urlopen(complete_url)
	myfile = f.read()
	return (myfile,unipro)

def get_location(protein_data):
	protein_list = protein_data[0].split('\n')[1:]
	protein_text =  ''.join(protein_list)
	motif_index = []
	for protein in range(0,len(protein_text)):
		motif_data =  protein_text[protein:protein+4]
		if (regex_check(motif_data) == 1):
			motif_index.append(str(protein + 1))

	return (protein_data[1], motif_index)
		

def regex_check(motif):
	# N{P}[ST]{P}
	# N and not P S or T not P
	my_regex_pattern = r'N[^P][S|T][^P]'
	if re.search(my_regex_pattern, motif):
		return 1
	else:
		return 0

if __name__ == '__main__':
	#dataset = open("/Users/krishna/MIRI/BSG/data/problem15.txt").read().split('\n')
	dataset = open("/Users/krishna/Downloads/rosalind_mprt.txt").read().split('\n')
	protein_data = []
	count = 0
	for unipro in dataset:
		protein_data.append(parse_url(unipro))

	for pd in protein_data:
		if (len(get_location(pd)[1]) > 0):
			print get_location(pd)[0]
			print ' '.join(get_location(pd)[1])