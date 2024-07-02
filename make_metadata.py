import pandas as pd
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description='Process fastq list file.')
parser.add_argument('-f', '--fastq_list', required=True, help='Path to the fastq list file')
args = parser.parse_args()

# Load the file
filelist = pd.read_csv(args.fastq_list)

# Set the column name
filelist.columns = ['sample']

# Extract only the relevant part of the sample names
filelist['sample'] = filelist['sample'].str.extract(r'raw_fastq/([\w\d_-]+_[\d_]+)\.fastq\.gz', expand=False)

# Extract drug and concentration information from the sample names
filelist['drug'] = filelist['sample'].str.split('_').str[1].astype(str)
filelist['conc_xMIC'] = filelist['sample'].str.split('_').str[2]

# Define the MIC labels mapping
mic_labels = {
    '03': 0.03125,
    '06': 0.0625,
    '12': 0.125,
    '25': 0.25,
    '37': 0.375,
    '50': 0.5,
    '75': 0.75,
    '1': 1,
    '2': 2
}

# Replace concentration labels with numerical values
filelist['conc_xMIC'] = filelist['conc_xMIC'].replace(mic_labels)

# Add a column for days_predepletion
filelist['days_predepletion'] = '1'

# Rename the columns
filelist.columns = ['column_name', 'drug', 'conc_xMIC', 'days_predepletion']

# Save the processed data to a new file
filelist.to_csv('sample_metadata.txt', sep='\t', index=False)

