#!/usr/bin/python3

import sys

def minimum_skew_positions (sequence):
	"""
	Minimum Skew Problem: Find a position in a genome minimizing the skew.
	Input: A DNA string Genome.
	Output: All integer(s) i minimizing Skew(Prefixi (Text)) among all values of i (from 0 to |Genome|).
	"""
	
	skew_vector = [0]
	skew_for_prefix = 0
	for i, nt in enumerate(sequence):
		if nt == 'C':
			skew_for_prefix -= 1
		elif nt == 'G':
			skew_for_prefix += 1
		skew_vector.append(skew_for_prefix)
	
	min_skew = min(skew_vector)
	minimum_skew_pos = [i for i, skew in enumerate(skew_vector) if skew == min_skew]
	
	return minimum_skew_pos
	
# Get filename from the command arguments and open the file
filename = str(sys.argv[1])
f = open(filename, 'r')

# Get the data from the file
sequence = f.readline().strip()

# Call the function that finds clumps
minimum_skew_pos = minimum_skew_positions(sequence)

print(' '.join(map(str,minimum_skew_pos)))
