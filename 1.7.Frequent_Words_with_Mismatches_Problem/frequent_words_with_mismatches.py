#!/usr/bin/python3

import sys

def frequent_words_with_mismatches(sequence, k, d):
	"""
	Frequent Words with Mismatches Problem: Find the most frequent k-mers with mismatches in a string.
	Input: A string Text as well as integers k and d. (You may assume k ≤ 12 and d ≤ 3.)
	Output: All most frequent k-mers with up to d mismatches in Text.
	"""
	
	kmer_counts = {}
	for i in range(len(sequence) - k + 1):
		kmer = sequence[i:i+k]
		if kmer not in kmer_counts:
			kmer_counts[kmer] = len(approximate_pattern_match_positions(kmer, sequence, d))
	
	max_count = max(kmer_counts.values())
	
	most_frequent_kmers = []
	for kmer, count in kmer_counts.items():
		if count == max_count:
			most_frequent_kmers.append(kmer)
	
	return most_frequent_kmers


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

# Get the data from the file
(sequence, pattern_width, allowed_mismatches) = f.readline().strip().split()

# Call function to calculate approximate pattern positions
most_frequent_kmers = frequent_words_with_mismatches(sequence, int(pattern_width), int(allowed_mismatches))
print(' '.join(most_frequent_kmers))