#!/usr/bin/python3

import sys
import operator

# Define a function to call later
def reverse_complement(seq):
	"""Return the reverse complement of a DNA string.""" 
	
	basecomplement = {'a':'t', 'c':'g', 't':'a', 'g':'c', 'A':'T', 'C':'G', 'T':'A', 'G':'C'}
	reverse_seq = seq[::-1]

	dna = '' # initialize the variable dna as an empty string
	for nt in reverse_seq:
		dna += basecomplement[nt] 
	return dna


# Get filename from the command arguments
filename = str(sys.argv[1])

# Open the file for reading
f = open(filename, 'r')

# Get the first line of the file - the DNA sequence
text = f.readline().strip()

# Call the funtion we defined earlier
revcom = reverse_complement(text)
print(revcom)