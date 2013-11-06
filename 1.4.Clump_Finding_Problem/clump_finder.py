#!/usr/bin/python3

import sys

def find_Kmers_forming_clumps (sequence, K, window, t):
	"""
	Clump Finding Problem: Find patterns forming clumps in a string.
	Input: A string Genome, and integers K, window, and t.
	Output: All distinct k-mers forming (window, t)-clumps in Genome.
	"""
	
	Kmers_forming_clumps = set() # Create a set. Duplicate Kmers are auto-discarded
	for i in range(len(sequence) - window + 1):
		window_seq = sequence[i:i+window]
		Kmers_forming_clumps.update(Kmers_found_more_than_t_times(sequence, K, t))

	return Kmers_forming_clumps
	
def Kmers_found_more_than_t_times(sequence, K, t):
	""" Find Kmers that are found at least t times in sequence """
	
	# Initialize a dictionary. Loop on the sequence and count Kmers
	counts = {}
	for i in range(len(sequence) - K + 1):
		Kmer = sequence[i:i+K]
		if Kmer not in counts:
			counts[Kmer] = 0
		counts[Kmer] += 1

	# Initialize a list to put accepted Kmers
	Kmers_more_than_t = []
	for Kmer, count in counts.items():
		if count >= t:
			Kmers_more_than_t.append(Kmer)

	return Kmers_more_than_t


# Get filename from the command arguments and open the file
filename = str(sys.argv[1])
f = open(filename, 'r')

# Get the data from the file
sequence = f.readline().strip()
(K, window, t) = map(int, f.readline().strip().split())

# Call the function that finds clumps
Kmers_forming_clumps = find_Kmers_forming_clumps(sequence, K, window, t)

# Print the results
print(' '.join(Kmers_forming_clumps))
