import re
from protein import *
import protein 

def code_protein(sequence):
	protein = []
	for i in xrange(0,len(sequence),3):
		if len(sequence[i:i+3]) > 2:
			protein.append(RNA_TABLE[sequence[i:i+3]])
	
	return ''.join(protein)

if __name__ == '__main__':
	data = 'AUGCUACUGGACGGACG'
	print data
	print code_protein(data)
