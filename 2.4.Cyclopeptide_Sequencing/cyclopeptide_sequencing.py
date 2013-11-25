#!/usr/bin/python3

import sys
import collections

class AminoAcidsMassTable:
	'''
	A class that serves as a wrapper for the amino acid mass table
	'''
	
	def __init__(self, mass_table_file):
		self.mass_dic = self._create_mass_dictionary_from_file(mass_table_file)
	
	def mass_for_amino_acid(self, aa):
		'''
		Return the mass for the requested amino acid
		'''
		return self.mass_dic[aa]
		
	def masses(self):
		'''
		Returns all unique masses for all amino acids
		'''
		return sorted(list(set(self.mass_dic.values())))
	
	def amino_acids(self):
		'''
		Returns all amino acids
		'''
		return sorted(self.mass_dic.keys())
	
	def _create_mass_dictionary_from_file(self, filename):
		'''
		Reads a tab delimited file with amino acid masses and returns 
		a dictionary with the masses.
		'''
		mass_dictionary = {}
		with open(filename) as f:
			for line in f:
				(aa, mass) = line.strip().split(" ")
				mass_dictionary[aa] = int(mass)
		return mass_dictionary

class Spectrum:
	'''
	A class that represents a mass spectrum.
	'''
	
	def __init__(self, masses):
		self.masses = sorted(masses)
	
	def is_identical_to(self, comparing_spectrum):
		'''
		Check if the spectrum is identical to another one.
		'''
		compare = lambda x, y: collections.Counter(x) == collections.Counter(y)
		#print (repr(self.masses) + "vs" + repr(comparing_spectrum.masses))
		return compare(self.masses, comparing_spectrum.masses)
	
	def is_consistent_with(self, comparing_spectrum):
		'''
		Check if the spectrum is consistent with a reference one.
		To be consistent, all its masses must be in the reference spectrum.
		'''
		for mass in self.masses:
			if mass not in comparing_spectrum.masses:
				return False
		return True
	
class Peptide:
	'''
	A class that represents a peptide.
	'''
	
	def __init__(self, aa_seq=None, aa_list=None, mass_list=None):
		self.aa_seq = aa_seq
		self.aa_list = aa_list
		self.mass_list = mass_list
		
		if self.aa_seq is not None and self.aa_list is None:
			self.aa_list = list(self.aa_seq)
		
	def theoretical_cyclospectrum(self, **args):
		'''
		Calculate the theoretical spectrum of a peptide given a table
		with the amino acid masses. Theoretical spectrum is the list
		of the total masses of all the cyclic subpeptides.
		'''
		if self.mass_list is None and 'aa_mass_table' not in args:
			raise Exception("FATAL")
		if 'aa_mass_table' in args:
			self.create_mass_list_from_aa_mass_table(aa_mass_table)
		
		theoretical_cyclospectrum = [0, self.total_mass()]
		for subpeptide in self.cyclic_subpeptides():
			theoretical_cyclospectrum.append(subpeptide.total_mass())
		
		return Spectrum(sorted(theoretical_cyclospectrum))
	
	def theoretical_linearspectrum(self, **args):
		'''
		Calculate the theoretical spectrum of a peptide given a table
		with the amino acid masses. Theoretical spectrum is the list
		of the total masses of all the cyclic subpeptides.
		'''
		if self.mass_list is None and 'aa_mass_table' not in args:
			raise Exception("FATAL")
		if 'aa_mass_table' in args:
			self.create_mass_list_from_aa_mass_table(aa_mass_table)
		
		theoretical_linearspectrum = [0, self.total_mass()]
		for subpeptide in self.subpeptides():
			theoretical_linearspectrum.append(subpeptide.total_mass())
		
		return Spectrum(sorted(theoretical_linearspectrum))
	
	def total_mass(self, **args):
		'''
		Calculate the total theoretical mass of a peptide given a table
		with the amino acid masses
		'''
		if self.mass_list is None and 'aa_mass_table' not in args:
			raise Exception("FATAL")
		if 'aa_mass_table' in args:
			self.create_mass_list_from_aa_mass_table(aa_mass_table)
		return sum(self.mass_list)
	
	def create_mass_list_from_aa_mass_table(self, aa_mass_table):
		'''
		Create a list with masses of the peptide, given a table with the
		amino acid masses
		'''
		if self.mass_list is None:
			self.mass_list = []
			for aa in self.aa_list:
				self.mass_list.append(aa_mass_table.mass_for_amino_acid(aa))
	
	def cyclic_subpeptides(self):
		'''
		Return all subpeptides for the Peptide.
		Also returns subpeptides that span the end and the start of the Peptide.
		'''
		for i in range(len(self)):
			for width in range(1, len(self)):
				subpeptide = self._extract_subpeptide_cyclically(width, i)
				yield subpeptide
	
	def subpeptides(self):
		'''
		Return all subpeptides for the Peptide.
		Also returns subpeptides that span the end and the start of the Peptide.
		'''
		for i in range(len(self)):
			for width in range(1, len(self)):
				subpeptide = self._extract_subpeptide(width, i)
				if subpeptide is not None:
					yield subpeptide
	
	def _extract_subpeptide(self, width, from_pos):
		'''
		Extract a part of the peptide from 'from_pos' with length 'width'. 
		If 'width' is larger than the peptide, continue from peptide start.
		'''
		subpeptide_aa_seq = ''
		subpeptide_aa_list = []
		subpeptide_mass_list = []
		for i in range(from_pos, from_pos+width):
			pos = i
			if pos >= len(self):
				return None
			if (self.aa_seq is not None):
				subpeptide_aa_seq += self.aa_seq[pos]
			if (self.aa_list is not None):
				subpeptide_aa_list.append(self.aa_list[pos])
			if (self.mass_list is not None):
				subpeptide_mass_list.append(self.mass_list[pos])
				
		return Peptide(aa_seq=subpeptide_aa_seq, aa_list=subpeptide_aa_list, mass_list=subpeptide_mass_list)
	
	def _extract_subpeptide_cyclically(self, width, from_pos):
		'''
		Extract a part of the peptide from 'from_pos' with length 'width'. 
		If 'width' is larger than the peptide, continue from peptide start.
		'''
		subpeptide_aa_seq = ''
		subpeptide_aa_list = []
		subpeptide_mass_list = []
		for i in range(from_pos, from_pos+width):
			pos = i
			if pos >= len(self):
				pos = pos - len(self)
			if (self.aa_seq is not None):
				subpeptide_aa_seq += self.aa_seq[pos]
			if (self.aa_list is not None):
				subpeptide_aa_list.append(self.aa_list[pos])
			if (self.mass_list is not None):
				subpeptide_mass_list.append(self.mass_list[pos])
				
		return Peptide(aa_seq=subpeptide_aa_seq, aa_list=subpeptide_aa_list, mass_list=subpeptide_mass_list)
	
	def __len__(self):
		if self.aa_seq is not None:
			return len(self.aa_seq)
		if self.aa_list is not None:
			return len(self.aa_list)
		if self.mass_list is not None:
			return len(self.mass_list)

def cyclopeptide_sequencing(spectrum, aa_mass_table):
	output_pep_list = []
	pep_list = [Peptide(mass_list=[])]
	while pep_list:
		pep_list = expand_list(pep_list, aa_mass_table.masses())
		index = 0
		while index >= 0 and index < len(pep_list):
			pep = pep_list[index]
			
			if (pep.theoretical_cyclospectrum().is_identical_to(spectrum)):
				output_pep_list.append(pep)
				pep_list.pop(index)
				continue
			
			if not pep.theoretical_linearspectrum().is_consistent_with(spectrum):
				#print ("Not consistent")
				#if (len(pep.theoretical_cyclospectrum().masses) < 3):
					#print (repr(pep.theoretical_cyclospectrum().masses))
				pep_list.pop(index)
				continue
			index += 1
	
	return output_pep_list
		
def expand_list(pep_list, aa_masses):
	'''
	Get a list of peptides (each peptide is defined as an array of its amino
	acid masses) and create a new list with the same peptides expanded
	by one amino acid.
	'''
	expanded_pep_list = []
	for pep in pep_list:
		for aa_mass in aa_masses:
			expanded_pep_list.append(Peptide(mass_list = pep.mass_list + [aa_mass]))
	return expanded_pep_list
	
# Get the command line arguments
dataset_file = sys.argv[1]
mass_table_file = sys.argv[2]

# Open the dataset file and read the spectrum masses
f = open(dataset_file)
input_masses = [int(mass) for mass in f.readline().strip().split(' ')]

# Create a spectrum
spectrum = Spectrum(input_masses)

# Create the amino acid mass table
aa_mass_table = AminoAcidsMassTable(mass_table_file)

# Run the algorithm
peptides_list = cyclopeptide_sequencing(spectrum, aa_mass_table)

# Print
mass_strings = []
for pep in peptides_list:
	mass_strings.append('-'.join(str(mass) for mass in pep.mass_list))
print (' '.join(mass_strings))

