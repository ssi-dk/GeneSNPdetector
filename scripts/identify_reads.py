#!/usr/bin/env python3

import os
import sys
import re

def get_read_paths_from_folder(directory):
	sample_read_dict = {}
	N_WGS_isolate_regex = "(?P<sample_name>.+?)(?P<sample_number>(_S[0-9]+)?)(?P<lane>(_L[0-9]+)?)_(?P<paired_read_number>R[1|2])(?P<set_number>(_[0-9]+)?)(?P<file_extension>\.fastq\.gz)"
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

sample_read_dict = get_read_paths_from_folder("/srv/data/MPV/MTG/S_epidermidis_Emeli/reads")
print(sample_read_dict)
		
