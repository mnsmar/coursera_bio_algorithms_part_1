#!/usr/bin/python3

import sys
import operator

# The filename is given as argument from the command line
# eg. > python frequent_words_counter.py filename1.txt
# The following command reads this argument
filename = str(sys.argv[1])

# Open the file for reading
f = open(filename, 'r')

# The first line is the text.  
# Note: strip removes the hidden character for newline from the end of the line
text = f.readline().strip()

# The second line is the Kmer size K
k = int(f.readline().strip())

# Initialize a dictionary. Loop on the text and count kmers
counts = {}
for i in range(len(text) - k + 1):
	kmer = text[i:i+k]
	if kmer not in counts:
		counts[kmer] = 0
	counts[kmer] += 1

# Find the maximum count
max_count = max(counts.values())

# Initialize a list to put the most frequent kmers
most_frequent_kmers = []
for kmer, count in counts.items():
	if count == max_count:
		most_frequent_kmers.append(kmer)

# Print the list with the most frequent kmers using space as separator
print(' '.join(most_frequent_kmers))