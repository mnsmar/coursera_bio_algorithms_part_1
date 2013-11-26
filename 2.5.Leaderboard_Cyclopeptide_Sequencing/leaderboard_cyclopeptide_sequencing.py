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
		self.stored_parent_mass = None
	
	def is_identical_to(self, comparing_spectrum):
		'''
		Check if the spectrum is identical to another one.
		'''
		compare = lambda x, y: collections.Counter(x) == collections.Counter(y)
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
	
	def parent_mass(self):
		'''
		Return the maximum mass as the parent mass of the spectrum
		'''
		if self.stored_parent_mass is None:
			self.stored_parent_mass = max(self.masses)
		return self.stored_parent_mass
	
class Peptide:
	'''
	A class that represents a peptide.
	'''
	
	def __init__(self, aa_seq=None, aa_list=None, mass_list=None):
		self.aa_seq = aa_seq
		self.aa_list = aa_list
		self.mass_list = mass_list
		self.stored_score = None
		self.stored_total_mass = None
		
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
		if self.stored_total_mass is None:
			if self.mass_list is None and 'aa_mass_table' not in args:
				raise Exception("FATAL")
			if 'aa_mass_table' in args:
				self.create_mass_list_from_aa_mass_table(aa_mass_table)
			self.stored_total_mass = sum(self.mass_list)
		
		return self.stored_total_mass
	
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
	
	def score(self, spectrum):
		'''
		Get the score of the peptide by comparing its spectrum to a reference one.
		The score is defined as the number of matching masses between the two spectra
		'''
		if self.stored_score is None:
			mass_dic = {}
			for mass in spectrum.masses:
				if mass in mass_dic:
					mass_dic[mass] += 1
				else:
					mass_dic[mass] = 1
			
			score = 0
			for mass in self.theoretical_cyclospectrum().masses:
				if mass in mass_dic and mass_dic[mass] > 0:
					score += 1
					mass_dic[mass] -= 1
		
			self.stored_score = score
			
		return self.stored_score
	
	def __len__(self):
		if self.aa_seq is not None:
			return len(self.aa_seq)
		if self.aa_list is not None:
			return len(self.aa_list)
		if self.mass_list is not None:
			return len(self.mass_list)

def leaderboard_cyclopeptide_sequencing(spectrum, aa_mass_table, N):
	leader_peptide = Peptide(mass_list=[])
	leaderboard = [leader_peptide]
	while leaderboard:
		leaderboard = expand_list(leaderboard, aa_mass_table.masses())
		index = 0
		while index >= 0 and index < len(leaderboard):
			pep = leaderboard[index]
			
			if pep.total_mass() == spectrum.parent_mass():
				if pep.score(spectrum) > leader_peptide.score(spectrum):
					leader_peptide = pep
			elif pep.total_mass() > spectrum.parent_mass():
				leaderboard.pop(index)
				continue
			index += 1
		leaderboard = cut(leaderboard, spectrum, N)
	
	print (leader_peptide.score(spectrum))
	return [leader_peptide]
		
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

def cut(leaderboard, spectrum, N):
	'''
	Return the top N highest scoring peptides including ties in leaderboard
	'''
	if len(leaderboard) < N:
		return leaderboard
	
	leaderboard.sort(key=lambda pep: pep.score(spectrum), reverse=True)
	
	min_score = leaderboard[N].score(spectrum)
	
	top_N_leaderboard = [pep for pep in leaderboard if pep.score(spectrum) >= min_score]
	return top_N_leaderboard
	
	
	
# Get the command line arguments
dataset_file = sys.argv[1]
mass_table_file = sys.argv[2]

# Open the dataset file and read the spectrum masses
f = open(dataset_file)
N = int(f.readline().strip())
input_masses = [int(mass) for mass in f.readline().strip().split(' ')]

#N = 26
#input_masses = [int(mass) for mass in "0 71 97 101 103 113 113 113 113 114 114 115 128 128 128 128 129 131 131 131 156 156 184 186 186 200 214 227 227 228 230 231 241 242 242 243 244 244 256 257 262 269 270 287 298 299 301 328 331 340 340 343 345 345 356 358 359 370 370 372 375 383 385 397 400 401 429 430 442 453 454 454 459 462 468 471 472 473 474 485 486 487 498 499 501 512 514 514 542 561 567 570 573 575 581 583 585 590 599 600 600 601 602 610 615 615 616 627 627 630 658 695 696 698 698 698 701 703 704 713 723 728 728 728 728 730 730 731 741 744 747 758 761 769 799 810 817 827 829 831 832 841 841 844 844 851 854 854 857 859 862 872 882 884 886 889 928 928 944 945 947 955 955 958 959 960 966 967 972 972 982 985 990 996 997 1000 1000 1003 1041 1056 1059 1062 1068 1068 1068 1073 1075 1075 1084 1087 1089 1095 1097 1103 1113 1114 1128 1128 1131 1152 1172 1172 1181 1182 1184 1189 1190 1190 1196 1197 1199 1200 1202 1210 1212 1227 1231 1242 1259 1259 1283 1295 1298 1303 1303 1303 1303 1304 1311 1312 1317 1318 1325 1325 1328 1330 1338 1340 1345 1355 1356 1388 1396 1416 1426 1426 1427 1431 1432 1432 1434 1440 1442 1443 1445 1451 1453 1453 1454 1458 1459 1459 1469 1489 1497 1529 1530 1540 1545 1547 1555 1557 1560 1560 1567 1568 1573 1574 1581 1582 1582 1582 1582 1587 1590 1602 1626 1626 1643 1654 1658 1673 1675 1683 1685 1686 1688 1689 1695 1695 1695 1696 1701 1703 1704 1713 1713 1733 1754 1757 1757 1771 1772 1782 1788 1790 1796 1798 1801 1810 1810 1812 1817 1817 1817 1823 1826 1829 1844 1882 1885 1885 1888 1889 1895 1900 1903 1913 1913 1918 1919 1925 1926 1927 1930 1930 1938 1940 1941 1957 1957 1996 1999 2001 2003 2013 2023 2026 2028 2031 2031 2034 2041 2041 2044 2044 2053 2054 2056 2058 2068 2075 2086 2116 2124 2127 2138 2141 2144 2154 2155 2155 2157 2157 2157 2157 2162 2172 2181 2182 2184 2187 2187 2187 2189 2190 2227 2255 2258 2258 2269 2270 2270 2275 2283 2284 2285 2285 2286 2295 2300 2302 2304 2310 2312 2315 2318 2324 2343 2371 2371 2373 2384 2386 2387 2398 2399 2400 2411 2412 2413 2414 2417 2423 2426 2431 2431 2432 2443 2455 2456 2484 2485 2488 2500 2502 2510 2513 2515 2515 2526 2527 2529 2540 2540 2542 2545 2545 2554 2557 2584 2586 2587 2598 2615 2616 2623 2628 2629 2641 2641 2642 2643 2643 2644 2654 2655 2657 2658 2658 2671 2685 2699 2699 2701 2729 2729 2754 2754 2754 2756 2757 2757 2757 2757 2770 2771 2771 2772 2772 2772 2772 2782 2784 2788 2814 2885".split(' ')]

# Create a spectrum
spectrum = Spectrum(input_masses)

# Create the amino acid mass table
aa_mass_table = AminoAcidsMassTable(mass_table_file)

# Run the algorithm
peptides_list = leaderboard_cyclopeptide_sequencing(spectrum, aa_mass_table, N)

# Print
mass_strings = []
for pep in peptides_list:
	mass_strings.append('-'.join(str(mass) for mass in pep.mass_list))
print (' '.join(mass_strings))

