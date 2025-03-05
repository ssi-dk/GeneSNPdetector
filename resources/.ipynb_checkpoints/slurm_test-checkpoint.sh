#!/bin/bash

#SBATCH -c 1                                # Number of CPU cores
#SBATCH --mem=4G                            # Memory allocation
#SBATCH -p outbreak                         # Partition name
#SBATCH --output=~/slurm_logs/sbatch-%j.log # Output log file (with job ID)
#SBATCH --job-name=debug_job
#SBATCH --output=debug_output.log
#SBATCH --error=debug_error.log

set -x  # Enable command tracing


# Load environment or other preparation steps
#source /users/home/thobec/.bashrc


# Command to execute
source /users/home/thobec/.bashrc; conda activate ariba; ariba run --force /dpssi/data/Projects/thej_Sepidermidis_point_resistance/proj/Tools/git.repositories/GeneSNPdetector/scripts/../resources/ariba_refs/S_epidermidis_AMR /srv/data/Projects/thej_Sepidermidis_point_resistance/proj/Resources/reads/RH_test/SH418x92_240603_NB501792_AHCCCWAFX7_R1.fastq.gz /srv/data/Projects/thej_Sepidermidis_point_resistance/proj/Resources/reads/RH_test/SH418x92_240603_NB501792_AHCCCWAFX7_R2.fastq.gz /dpssi/data/Projects/thej_Sepidermidis_point_resistance/proj/Analysis/ariba/RH_test/SH418x92_240603_NB501792_AHCCCWAFX7