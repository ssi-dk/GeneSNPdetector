#!/usr/bin/env python3

import os
import sys
import argparse

import SNPfinder as sf


def cmdline_args():
	p = argparse.ArgumentParser(description="Idenity specified SNPs in all queried sequences")
	p.add_argument("Query", nargs="*",
					help = "Contig files to query")
	p.add_argument("-o", "--output", default = "SNPfinder_out",
					help="Output directory")
	p.add_argument("-s", "--scheme",
					help="Name of scheme to check")
	p.add_argument("-m", "--matrix_output", action="store_true", default=True,
					help ="Output presence/absence matrix of SNPs in reference")
	p.add_argument("-l", "--list_output", action='store_true', default=True,
					help ="Output list of reference SNPs found")

	return(p.parse_args())



def main(args):
	snps_dict = {}
	snp_lists_dict = {}
	if not os.path.exists(args.output):
		os.makedirs(args.output)
	SNPtable = os.path.join(os.path.dirname(os.path.realpath(__file__)),"..","resources",args.scheme+".tsv")
	reference_fasta = os.path.join(os.path.dirname(os.path.realpath(__file__)),"..","resources",args.scheme+".fasta")
	snp_reference,lineage_dict,lineage_list = sf.load_reference_snps(SNPtable)

	for query in args.Query:
		ID = ".".join(query.split('/')[-1].split('.')[:-1])
		prefix = os.path.join(args.output,ID)
		sf.run_nucmer_and_showsnps(query,reference_fasta,prefix)
		snps,snp_list = sf.extract_from_snpfile(prefix,snp_reference)
		snps_dict[ID] = snps
		snp_lists_dict[ID] = snp_list
	matrix_out_file = os.path.join(args.output,"SNPfinder_binary_matrix.tsv")
	list_out_file = os.path.join(args.output,"SNPfinder_list.tsv")
	lineage_out_file = os.path.join(args.output,"SNPfinder_type_SNPs.txt")
	#sf.print_matrix(snps_dict,snp_reference,matrix_out_file)
	sf.print_list(snp_lists_dict,snp_reference,list_out_file)
	sf.print_matrix_and_type(snps_dict,snp_reference,matrix_out_file,lineage_dict,lineage_list,lineage_out_file)

		


if __name__ == '__main__':
	args = cmdline_args()
	main(args)