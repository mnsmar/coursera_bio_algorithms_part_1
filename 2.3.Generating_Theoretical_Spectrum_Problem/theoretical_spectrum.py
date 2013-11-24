#!/usr/bin/python3

import sys

def create_theoretical_spectrum(peptide, mass_table):
	theoretical_spectrum = [0, peptide_mass(peptide, mass_table)]
	for i in range(len(peptide)):
		for width in range(1,len(peptide)):
			fragment = extract_fragment_from_cyclic_seq(peptide, width, i)
			fragment_mass = peptide_mass(fragment, mass_table)
			theoretical_spectrum.append(fragment_mass)
	theoretical_spectrum.sort()
	return theoretical_spectrum

def extract_fragment_from_cyclic_seq(seq, fragment_width, from_pos):
	'''
	Given a linear seq representing a cyclic seq, extract a sequence fragment
	of given length starting at a given position
	'''
	fragment = ''
	for i in range(from_pos, from_pos+fragment_width):
		if i < len(seq):
			fragment += seq[i]
		else:
			fragment += seq[i-len(seq)]
	return fragment

def peptide_mass(peptide, mass_table):
	'''
	Calculate the theoretical mass of a peptide given a table
	with the amino acid masses
	'''
	mass = 0
	for aa in peptide:
		mass += mass_table[aa]
	return mass

def create_mass_dictionary_from_file(filename):
	mass_dictionary = {}
	with open(filename) as f:
		for line in f:
			(aa, mass) = line.strip().split(" ")
			mass_dictionary[aa] = int(mass)
	return mass_dictionary

dataset_file = sys.argv[1]
mass_table_file = sys.argv[2]

f = open(dataset_file)
peptide = f.readline().strip()

mass_dictionary = create_mass_dictionary_from_file(mass_table_file)
theoretical_spectrum = create_theoretical_spectrum(peptide, mass_dictionary)

print (' '.join([str(i) for i in theoretical_spectrum]))


