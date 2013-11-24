#!/usr/bin/python3

import sys

def create_codon_dictionary_from_file(filename):
	codon_dictionary = {}
	with open(filename) as f:
		for line in f:
			if " \n" in line:
				triplet = line.strip()
				codon_dictionary[triplet] = "stop"
				triplet = triplet.replace("U", "T")
				codon_dictionary[triplet] = "stop"
			else:
				(triplet, aa) = line.strip().split(" ")
				codon_dictionary[triplet] = aa
				triplet = triplet.replace("U", "T")
				codon_dictionary[triplet] = aa
	return codon_dictionary

def translate_to_protein(rna, codon_dictionary):
	protein = ''
	for i in range(0, len(rna)-2, 3):
		codon = rna[i:i+3]
		if codon_dictionary[codon] != "stop":
			protein += codon_dictionary[codon]
		else:
			break
	return protein

def reverse_complement(seq):
	"""Return the reverse complement of a DNA string.""" 
	
	basecomplement = {'a':'t', 'c':'g', 't':'a', 'g':'c', 'u':'a', 'A':'T', 'C':'G', 'T':'A', 'G':'C', 'U':'A'}
	reverse_seq = seq[::-1]

	dna = '' # initialize the variable dna as an empty string
	for nt in reverse_seq:
		dna += basecomplement[nt] 
	return dna

def patterns_encoding_for_peptide_in_rna(rna, peptide):
	patterns_encoding_for_peptide = []
	for i in range(len(rna) - 3*len(peptide) + 1):
		rna_part = rna[i:i+3*len(peptide)]
		test_peptide = translate_to_protein(rna_part, codon_dictionary)
		if test_peptide == peptide:
			patterns_encoding_for_peptide.append(rna_part)
	return patterns_encoding_for_peptide


dataset_file = sys.argv[1]
codon_table_file = sys.argv[2]

f = open(dataset_file)
dna = f.readline().strip()
peptide = f.readline().strip()

rev_comp_dna = reverse_complement(dna)
codon_dictionary = create_codon_dictionary_from_file(codon_table_file)

patterns_encoding_for_peptide = patterns_encoding_for_peptide_in_rna(dna, peptide)
for pattern in patterns_encoding_for_peptide_in_rna(rev_comp_dna, peptide):
	patterns_encoding_for_peptide.append(reverse_complement(pattern))
print ("\n".join(patterns_encoding_for_peptide))



