#!/usr/bin/python3

import sys

def approximate_pattern_match_positions(pattern, sequence, d):
	"""
	Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.
	   Input: Two strings Pattern and Text along with an integer d.
	   Output: All positions where Pattern appears in Text with at most d mismatches.
	"""
	
	pattern_length = len(pattern)
	pattern_positions = []
	for i in range(len(sequence) - pattern_length + 1):
		kmer = sequence[i:i+pattern_length]
		if edit_distance(pattern, kmer) <= d:
			pattern_positions.append(str(i))
			
	return pattern_positions


def edit_distance(pattern1, pattern2):
	"""
	Compare two patterns and calculate edit distance - the number of mismatches between the sequences
	"""
	
	edit_distance = 0
	for i, nt in enumerate(pattern1):
		if (nt != pattern2[i]):
			edit_distance += 1
	
	return edit_distance

# Get filename from the command arguments and open the file
filename = str(sys.argv[1])
f = open(filename, 'r')

# The first line in the file is the pattern
pattern = f.readline().strip()
sequence = f.readline().strip()
allowed_mismatches = int(f.readline().strip())

# Call function to calculate approximate pattern positions
pattern_positions = approximate_pattern_match_positions(pattern, sequence, allowed_mismatches)
print(' '.join(pattern_positions))