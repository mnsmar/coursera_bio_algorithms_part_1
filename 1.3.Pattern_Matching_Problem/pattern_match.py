#!/usr/bin/python3

import sys
import operator

# Get filename from the command arguments
filename = str(sys.argv[1])

# Open the file for reading
f = open(filename, 'r')

# The first line in the file is the pattern
pattern = f.readline().strip()
pattern_length = len(pattern)

# The second line is the file is the sequence
sequence = f.readline().strip()

# Loop on the sequence and check if pattern is found.
# Before the loop initialize an array to store the positions found.
pattern_positions = []
for i in range(len(sequence) - pattern_length + 1):
	kmer = sequence[i:i+pattern_length]
	if kmer == pattern:
		pattern_positions.append(str(i))

# Print the list with the positions
print(' '.join(pattern_positions))