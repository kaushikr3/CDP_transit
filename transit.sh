#!/bin/bash -l
#SBATCH --partition=scu-cpu   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=cgi_42
#SBATCH --time=24:00:00   # HH/MM/SS
#SBATCH --mem=8G   # memory requested, units available: K,M,G,T
#SBATCH --output=job_output_%j.log
#SBATCH --error=job_error_%j.log

# Load the conda env
source ~/.bashrc
conda activate transit

echo '%%%%% ND base RUN %%%%%'
# run the make metadata script to make the metadata sheet for the CDP run
python ~/CDP/scripts/make_metadata.py -f fastq_list.txt

# Create temp_raw_fastq directory for transit to run
mkdir temp_raw_fastq

## RUN the run_transit.sh script

# -i is the list of fastq files in the run(has to be generated)
# -r is the sgrna info file, doesn't change unless library changes
# -m is the sample metadata sheet that has to be generated
# -o is the output prefix(text input)
# -f is the flag to run the BSAF comparisons
# -b is the base name for the BSAF run like ND_bsaf, NDBSAF, BSAF etc.

sh ~/CDP/scripts/run_transit.sh -i fastq_list.txt -r ~/CDP/sgRNA_info_RLC0025.txt -m sample_metadata.txt -o CDP42

# BSAF RUN:

echo '%%%%% bsAf base RUN %%%%%'
# run the make metadata script to make the metadata sheet for the CDP run
python ~/CDP/scripts/make_metadata.py -f fastq_list_bsaf.txt

## RUN the run_transit.sh script
sh ~/CDP/scripts/run_transit.sh -i fastq_list.txt -r ~/CDP/sgRNA_info_RLC0025.txt -m sample_metadata.txt -o CDP42 -f -b NDBSAF
