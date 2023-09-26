#!/usr/bin/env python3

import os
import sys


def run_nucmer(query,reference,prefix):
	cmd = "nucmer --prefix={prefix} {reference} {query}".format(prefix=prefix,reference=reference,query=query)
	os.system(cmd)

def run_showsnps(delta,snps_out):
	cmd = "show-snps -HrT {} > {}".format(delta,snps_out)
	os.system(cmd)

def run_nucmer_and_showsnps(query,reference,prefix):
	cmd = "nucmer --prefix={prefix} {reference} {query}; show-snps -HrT {prefix}.delta > {prefix}.snps".format(prefix=prefix,reference=reference,query=query)
	os.system(cmd)

def load_reference_snps(reference_snp_file):
	snp_reference = {}
	type_dict = {}
	type_list = []
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
					type_idx = line.index("type")
				except:
					print("Incorrect headers in snp reference file\tMust contain \"ref_name\", \"position\", \"ref_base\", \"alt_base\"")
			else:
				try: 
					snp_reference[line[ref_idx]][line[pos_idx]] = [line[refbase_idx],line[altbase_idx],line]
				except:
					snp_reference[line[ref_idx]] = {line[pos_idx]:[line[refbase_idx],line[altbase_idx],line]}
				type_dict[line[ref_idx]+"::"+line[refbase_idx]+line[pos_idx]+line[altbase_idx]] = line[type_idx]
				if not line[type_idx] in type_list:
					type_list.append(line[type_idx])
	return(snp_reference,type_dict,type_list)


def extract_all_from_snpfile(prefix):
	showsnps_file = prefix+".snps"
	snps = []
	snp_list = []
	with open(showsnps_file) as f:
		for line in f:
			line = line.rstrip('\n').split('\t')
			ref_name = line[10]
			pos = line[0]
			ref_base = line[1]
			alt_base = line[2]
			snps.append(ref_name+"::"+ref_base+pos+alt_base)
			snp_list.append(snp_reference[ref_name][pos][2])
	return(snps,snp_list)

def extract_from_snpfile(prefix,snp_reference):
	showsnps_file = prefix+".snps"
	snps = []
	snp_list = []
	with open(showsnps_file) as f:
		for line in f:
			line = line.rstrip('\n').split('\t')
			ref_name = line[10]
			if ref_name in snp_reference:
				pos = line[0]
				if pos in snp_reference[ref_name]:
					ref_base = line[1]
					alt_base = line[2]
					query_name = line[11]
					query_pos = line[3]
					snps.append(ref_name+"::"+ref_base+pos+alt_base)
					snp_list.append([ref_name]+[ref_name+"::"+ref_base+pos+alt_base]+snp_reference[ref_name][pos][2][1:]+[query_name,query_pos])
					print(snp_list)
	return(snps,snp_list)


def print_matrix(snps_dict,snp_reference,matrix_out_file):
	header = []
	for ref_name in snp_reference:
		for pos in snp_reference[ref_name]:
			snp = ref_name+"::"+snp_reference[ref_name][pos][0]+pos+snp_reference[ref_name][pos][1]
			header.append(snp)
	o = open(matrix_out_file,"w")
	o.write("query\t"+"\t".join(header)+"\n")
	for ID in snps_dict:
		printline = [ID]
		for snp in header:
			if snp in snps_dict[ID]:
				printline.append("1")
			else:
				printline.append("0")
		o.write("\t".join(printline)+"\n")


def print_list(snp_lists_dict,snp_reference,list_out_file):
	o = open(list_out_file,"w")
	o.write("query\tref_name\tSNP\tpos\tref_base\talt_base\tnotes\tDOI\tarticle_link\tquery_name\tquery_pos\n")
	for ID in snp_lists_dict:
		for snp_list in snp_lists_dict[ID]:
			print("\t".join(snp_list))
			o.write(ID+"\t"+"\t".join(snp_list)+"\n")



def print_matrix_and_type(snps_dict,snp_reference,matrix_out_file,type_dict,type_list,type_out_file):
	header = []
	for ref_name in snp_reference:
		for pos in snp_reference[ref_name]:
			snp = ref_name+"::"+snp_reference[ref_name][pos][0]+pos+snp_reference[ref_name][pos][1]
			header.append(snp)
	o = open(matrix_out_file,"w")
	o.write("query\t"+"\t".join(header)+"\n")

	o2 = open(type_out_file,"w")
	o2.write("query\t"+"\t".join(type_list)+"\n")

	for ID in snps_dict:
		printline = [ID]
		type_counts = {}
		for type in type_list:
			type_counts[type] = 0
		for snp in header:
			if snp in snps_dict[ID]:
				printline.append("1")
				type_counts[type_dict[snp]] += 1
			else:
				printline.append("0")
		o.write("\t".join(printline)+"\n")
		o2_line = [ID]
		for type in type_list:
			o2_line.append(str(type_counts[type]))
		o2.write("\t".join(o2_line)+"\n")
	o.close()
	o2.close()