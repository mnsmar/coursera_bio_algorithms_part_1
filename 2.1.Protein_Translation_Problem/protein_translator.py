# Protein Translation Problem: Translate an RNA string into an amino acid string.
# Input: An RNA string Pattern.
# Output: The translation of Pattern into an amino acid string Peptide.

# python ex1.py ex1.txt RNA_codon_table_1.txt

import sys

rna_file = sys.argv[1]
codon_table_file = sys.argv[2]

rna_txt = open(rna_file)
rna_seq = rna_txt.readline().strip()

codon_dictionary = {}
with open(codon_table_file) as f:
	for line in f:
		if " \n" in line:
			triplet = line.strip()
			codon_dictionary[triplet] = "stop"
		else:
			(triplet, aa) = line.strip().split(" ")
			codon_dictionary[triplet] = aa

protein = ''
for i in xrange(0,len(rna_seq)-2,3):
	codon = rna_seq[i:i+3]
	if codon_dictionary[codon] != "stop":
		protein += codon_dictionary[codon]
	else:
		break

print protein
