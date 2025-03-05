#!/usr/bin/env python3

import os
import sys


def run_nucmer(query,reference,prefix):
	cmd = "nucmer --prefix={prefix} {reference} {query}".format(prefix=prefix,reference=reference,query=query)
	print(cmd)
	#os.system(cmd)

def run_showsnps(delta,snps_out):
	cmd = "show-snps -HrT {} > {}".format(delta,snps_out)
	print(cmd)
	#os.system(cmd)

def run_nucmer_and_showsnps(query,reference,prefix):
	cmd = "nucmer --prefix={prefix} {reference} {query}; show-snps -HrT {prefix}.delta > {prefix}.snps".format(prefix=prefix,reference=reference,query=query)
	print(cmd)
	os.system(cmd)

def load_reference_snps(reference_snp_file):
	snp_reference = {}
	with open(reference_snp_file) as f:
		firstline = True
		for line in f:
			line = line.strip('\n').split('\t')
			if firstline:
				firstline = False
				try:
					ref_idx = line.index("ref_name")
					pos_idx = line.index("position")
					refbase_idx = line.index("ref_base")
					altbase_idx = line.index("alt_base")
				except:
					print("Incorrect headers in snp reference file\tMust contain \"ref_name\", \"position\", \"ref_base\", \"alt_base\"")
			else:
				try: 
					snp_reference[line[ref_idx]][line[pos_idx]] = [line[refbase_idx],line[altbase_idx]]
				except:
					snp_reference[line[ref_idx]] = {line[pos_idx]:[line[refbase_idx],line[altbase_idx]]}
	return(snp_reference)


def check_showsnps(prefix,snp_reference):
	showsnps_file = prefix+".snps"
	snps = []
	with open(showsnps_file) as f:
		for line in f:
			line = line.rstrip('\n').split('\t')
			ref_name = line[10]
			print(ref_name)
			if ref_name in snp_reference:
				pos = line[0]
				print(pos)
				if pos in snp_reference[ref_name]:
					ref_base = line[1]
					alt_base = line[2]
					snps.append(ref_name+"::"+ref_base+pos+alt_base)
	return(snps)




run_nucmer_and_showsnps("../assemblies/GAS-2023-0003.fasta","../studies/genes/virulence/speA_speC_ssrA_rofA.fasta","ref_qryUK")

snp_reference = load_reference_snps("/srv/data/tools/git.repositories/GeneSNPdetector/resources/SNP_reference_table.tsv")
print(snp_reference)
snps = check_showsnps("ref_qryUK",snp_reference)
print(snps)