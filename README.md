# CDP_transit
Chemical drug profiling pipeline using TRANSIT2

			
## run_transit.sh
Bash script that runs the TRANSIT2 pipeline for CDP data (FASTQ files).

## transit.sh
The slurm script to run the pipeline. 

Required arguments:

-i is the list of fastq files in the run. 

-r is the sgrna info file for the library of interest.

-m is the sample metadata sheet that has to be generated.

-o is the output prefix(text input).

-f is the flag to run the BSAF comparisons.

-b is the base name for the BSAF run like ND_bsaf, NDBSAF, BSAF etc. (only required with -f).

## extract_drugs.py
Python script that extracts the names of the Drugs used in the experiment from the FASTQ files.

Required arguments:

sample_metadata file that holds the information about the drugs and their concentrations.

## make_metadata.py
Python script to make the sample_metadata file that is required for TRANSIT2.

Required arguments:

-f, --fastq_list: text file with the list of fastq files in the experiment.

## add_rvnumber.py
Python script that edits the results file generated by TRANSIT2 to add Rv_number column.

Required arguments:

-d,--data: The results file generated from TRANSIT in tsv format.

-o,--output: The name of the output file.

run_transit_ND.sh,run_transit_nd_only.sh,run_transit_bsaf_only.sh,run_transit_both.sh: scripts that were created while testing
