#!/usr/bin/env python3

import os
import sys
import re
import subprocess
import argparse
my_env = os.environ


def parse_args(argv):
	parser = argparse.ArgumentParser(description='Run LEG typing (SBT typing) on LEG samples from a given year')
	parser.add_argument("-i", "--input",
						help = "Input directory containing read files with standard paired end Illumina naming",
						type=str,
						required = True)
	parser.add_argument("-o", "--output",
						help = "Output directory. Will contain a folder with ariba output for each file as well as a summary text file.",
						type=str,
						required = True)
	parser.add_argument("-r", "--reference",
						help = "Name of ariba reference located in the resources folder of the github directory.",
						type=str,
						required = True)
	parser.add_argument("-p", "--partition",
						help = "Name of ariba reference located in the resources folder of the github directory.",
						type=str,
						default = "project")
	args = parser.parse_args()
	return args

def get_read_paths_from_folder(directory):
	sample_read_dict = {}
	#N_WGS_isolate_regex = "(?P<sample_name>.+?)(?P<sample_number>(_S[0-9]+)?)(?P<lane>(_L[0-9]+)?)_(?P<paired_read_number>R[1|2])(?P<set_number>(_[0-9]+)?)(?P<file_extension>\.fastq\.gz)"
	N_WGS_isolate_regex = "(?P<sample_name>.+?)\.(?P<paired_read_number>R[1|2])(?P<file_extension>\.fastq\.gz)"
	read_regex = re.compile(N_WGS_isolate_regex)
	files = os.listdir(directory)
	for file in files:
		read_match = re.match(read_regex,file)
		if read_match:
			sample_name = read_match.group("sample_name")
			if not sample_name in sample_read_dict:
				if read_match.group("paired_read_number") == "R1":
					sample_read_dict[sample_name] = [os.path.join(directory,file),""]
				elif read_match.group("paired_read_number") == "R2":
					sample_read_dict[sample_name] = ["",os.path.join(directory,file)]
			else:
				if read_match.group("paired_read_number") == "R1":
					print("hey")
					if sample_read_dict[sample_name][0] == "":
						sample_read_dict[sample_name][0] = os.path.join(directory,file)
				elif read_match.group("paired_read_number") == "R2":
					if sample_read_dict[sample_name][1] == "":
						sample_read_dict[sample_name][1] = os.path.join(directory,file)
	return(sample_read_dict)

def queue_ariba(ariba_reference,r1_path,r2_path,output_dir,partition):
	cmd = "ariba run --force {ariba_reference} {r1_path} {r2_path} {output_dir}".format(ariba_reference=ariba_reference,r1_path=r1_path,r2_path=r2_path,output_dir=output_dir)
	print(cmd)
	os.system(cmd)
	#sbatch_cmd = "sbatch -D . -c 4 --mem=12G -J ariba_var_caller -p {partition} -wrap=\"{cmd}\"".format(partition=partition,cmd=cmd)
	#print(sbatch_cmd)
	#process = subprocess.Popen("sbatch -D . -c 4 --mem=12G -J ariba_var_caller -p {partition} --wrap=\'{cmd}\'".format(partition=partition,cmd=cmd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, env=my_env, encoding='utf-8')
	#process_out, process_err = process.communicate()
	#sys.stdout.write("{}\n".format(process_out))
	#sys.stderr.write("{}\n".format(process_err))
	#job_id = process_out.split(' ')[-1].strip()
	#return(job_id)
	return(cmd)

def summarize_ariba_results(ariba_report_path,ariba_reference,output_file,partition,dependency_job_ids):
	cmd = "source activate env_r-markdown; Rscript /srv/data/tools/git.repositories/GeneSNPdetector/scripts/ariba_summary.R {ariba_report_path} {output_file} {ariba_reference}; conda deactivate".format(ariba_report_path=ariba_report_path,output_file=output_file,ariba_reference=ariba_reference)
	print(cmd)
	if len(dependency_job_ids) == 0:
		process = subprocess.Popen("sbatch -D . -c 1 --mem=4G -J ariba_var_caller_summary -p {partition} --wrap=\'{cmd}\'".format(partition=partition,cmd=cmd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, env=my_env, encoding='utf-8')
	else:
		process = subprocess.Popen("sbatch -D . -c 1 --mem=4G -J ariba_var_caller_summary -p {partition} --dependency=afterany:{dependency_job_ids} --wrap=\'{cmd}\'".format(partition=partition,dependency_job_ids=dependency_job_ids,cmd=cmd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, env=my_env, encoding='utf-8')
	process_out, process_err = process.communicate()
	sys.stdout.write("{}\n".format(process_out))
	sys.stderr.write("{}\n".format(process_err))
	job_id = process_out.split(' ')[-1].strip()
	return(job_id)


args = parse_args(sys.argv)

sample_read_dict = get_read_paths_from_folder(args.input)
print(sample_read_dict)

if not os.path.exists(args.output):
	os.makedirs(args.output)

reference_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),"..","resources","ariba_refs",args.reference)
for sample_name in sample_read_dict:
	job_ids = []
	sample_ariba_out_dir = os.path.join(args.output,sample_name)
	if not os.path.exists(os.path.join(sample_ariba_out_dir,"report.tsv")):
		job_id = queue_ariba(reference_dir,sample_read_dict[sample_name][0],sample_read_dict[sample_name][1],sample_ariba_out_dir,args.partition)
		job_ids.append(job_id)
summary_job_id = summarize_ariba_results(args.output,args.reference,os.path.join(args.output,"ariba_summary.tsv"),args.partition,job_ids)
print(job_ids)

