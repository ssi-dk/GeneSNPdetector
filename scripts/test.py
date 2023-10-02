#!/usr/bin/env python3
import os
import sys
from math import floor

def get_amino_acid_mutation(nt_seq,nt_pos,nt_ref,nt_alt):
	#DNA_Nucleotides = ['A', 'C', 'G', 'T']
	#DNA_ReverseComplement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
	DNA_Codons = {
		# 'M' - START, '_' - STOP
		"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
		"TGT": "C", "TGC": "C",
		"GAT": "D", "GAC": "D",
		"GAA": "E", "GAG": "E",
		"TTT": "F", "TTC": "F",
		"GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
		"CAT": "H", "CAC": "H",
		"ATA": "I", "ATT": "I", "ATC": "I",
		"AAA": "K", "AAG": "K",
		"TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
		"ATG": "M",
		"AAT": "N", "AAC": "N",
		"CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
		"CAA": "Q", "CAG": "Q",
		"CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
		"TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
		"ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
		"GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
		"TGG": "W",
		"TAT": "Y", "TAC": "Y",
		"TAA": "_", "TAG": "_", "TGA": "_"
	}
	codon_start_idx = floor((nt_pos-1)/3)*3
	aa_pos = floor((nt_pos-1)/3)+1
	codon_mut_idx = nt_pos-codon_start_idx-1
	ref_codon = nt_seq[codon_start_idx:(codon_start_idx+3)]
	alt_codon = list(ref_codon)
	alt_codon[codon_mut_idx] = nt_alt
	alt_codon = "".join(alt_codon)
	print(ref_codon)
	aa_ref = DNA_Codons[ref_codon]
	aa_alt = DNA_Codons[alt_codon]
	return(aa_pos,aa_ref,aa_alt)


aa_pos,aa_ref,aa_alt = get_amino_acid_mutation("ATGGCTTGC",5,"C","A")
print(aa_pos)
print(aa_ref)
print(aa_alt)