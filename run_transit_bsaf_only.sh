#!/bin/bash

# Usage: ./run_transit_bsaf.sh -i input_filelist.txt -r sgRNA_info.txt -m sample_metadata.txt -o output_prefix

# Load necessary modules or environment
spack load /kinvl2y # r@4.2.2%gcc@8.2.0

# Parse arguments
while getopts "i:r:m:o:" opt; do
  case $opt in
    i) input_filelist="$OPTARG"
    ;;
    r) sgRNA_info="$OPTARG"
    ;;
    m) sample_metadata="$OPTARG"
    ;;
    o) output_prefix="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

# Check if input filelist is provided
if [ -z "$input_filelist" ]; then
  echo "Input filelist is required. Usage: ./run_transit_bsaf.sh -i input_filelist.txt -r sgRNA_info.txt -m sample_metadata.txt -o output_prefix"
  exit 1
fi

# Check if sgRNA info file is provided
if [ -z "$sgRNA_info" ]; then
  echo "sgRNA info file is required. Usage: ./run_transit_bsaf.sh -i input_filelist.txt -r sgRNA_info.txt -m sample_metadata.txt -o output_prefix"
  exit 1
fi

# Check if sample metadata is provided
if [ -z "$sample_metadata" ]; then
  echo "Sample metadata is required. Usage: ./run_transit_bsaf.sh -i input_filelist.txt -r sgRNA_info.txt -m sample_metadata.txt -o output_prefix"
  exit 1
fi

# Check if output prefix is provided
if [ -z "$output_prefix" ]; then
  echo "Output prefix is required. Usage: ./run_transit_bsaf.sh -i input_filelist.txt -r sgRNA_info.txt -m sample_metadata.txt -o output_prefix"
  exit 1
fi

# Create the counts_files directory
mkdir -p counts_files
counts_dir="counts_files"

# Function to extract counts for each FASTQ file
extract_counts() {
  fastq_file=$1
  base_name=$(basename "$fastq_file" .fastq.gz)
  output_file="${counts_dir}/${base_name}.counts"
  echo "Extracting counts for $fastq_file"
  transit cgi extract_counts "$fastq_file" "$sgRNA_info" "Barcode" "$output_file" -delete_temp_fastQ
}

# Read the input file list and extract counts for each file
count_files=()
while IFS= read -r fastq_file; do
  extract_counts "$fastq_file"
  base_name=$(basename "$fastq_file" .fastq.gz)
  count_file="${counts_dir}/${base_name}.counts"
  count_files+=("$count_file")
done < "$input_filelist"

# Extract unique drug names from sample metadata
unique_drugs=$(python extract_drugs.py "$sample_metadata")

# Create combined counts files for ND_bsAf base only
nd_bsaf_files=($(printf "%s\n" "${count_files[@]}" | grep "ND_bsAf" | sort))

nd_bsaf_headers=()
for file in "${nd_bsaf_files[@]}"; do
  nd_bsaf_headers+=($(basename "$file" .counts))
done

# Create combined counts file for ND_bsAf
if [ ${#nd_bsaf_files[@]} -gt 0 ]; then
  combined_counts_nd_bsaf="${output_prefix}_ND_bsaf_combined_counts.txt"
  echo "Creating combined counts for ND_bsAf base"
  header_str_nd_bsaf=$(IFS=,; echo "${nd_bsaf_headers[*]}")
  count_files_str_nd_bsaf="${nd_bsaf_files[*]}"
  transit cgi create_combined_counts "$header_str_nd_bsaf" $count_files_str_nd_bsaf "$combined_counts_nd_bsaf"
fi

# Function to generate the create_combined_counts command for each drug
create_combined_counts_for_drug() {
  drug=$1
  drug_files=($(printf "%s\n" "${count_files[@]}" | grep "$drug" | sort))
  if [ ${#drug_files[@]} -gt 0 ]; then
    headers=()
    for file in "${drug_files[@]}"; do
      headers+=($(basename "$file" .counts))
    done
    
    # Create combined headers and files for ND_bsAf base
    headers_nd_bsaf=("${headers[@]}" "${nd_bsaf_headers[@]}")
    combined_files_nd_bsaf=("${drug_files[@]}" "${nd_bsaf_files[@]}")
    combined_counts_nd_bsaf="${output_prefix}_${drug}_combined_counts_nd_bsaf.txt"
    if [ ${#nd_bsaf_files[@]} -gt 0 ]; then
      echo "Creating combined counts for $drug with ND_bsAf base"
      header_str_nd_bsaf=$(IFS=,; echo "${headers_nd_bsaf[*]}")
      count_files_str_nd_bsaf="${combined_files_nd_bsaf[*]}"
      transit cgi create_combined_counts "$header_str_nd_bsaf" $count_files_str_nd_bsaf "$combined_counts_nd_bsaf"
    fi
  fi
}

# Generate create_combined_counts command for each drug
for drug in $unique_drugs; do
  create_combined_counts_for_drug "$drug"
done

# Extract fractional abundances for each combined counts file
fractional_abundances_nd_bsaf="${output_prefix}_frac_abundances_nd_bsaf.txt"

for drug in $unique_drugs; do
  combined_counts_nd_bsaf="${output_prefix}_${drug}_combined_counts_nd_bsaf.txt"
  
  echo "Extracting fractional abundances for $drug with ND_bsAf base"
  transit cgi extract_abund "$combined_counts_nd_bsaf" "$sample_metadata" ND_bsAf "$sgRNA_info" Pred_logFC ORF "$drug" 1 "${output_prefix}_${drug}_frac_abundances_nd_bsaf.txt" -no_uninduced #change the ND_bsaf argument based on the bsaf name of the files.
done

# Run CRISPRi-DR model for each fractional abundance file
for drug in $unique_drugs; do
  fractional_abundances_nd_bsaf="${output_prefix}_${drug}_frac_abundances_nd_bsaf.txt"
  
  crispri_results_nd_bsaf="${output_prefix}_${drug}_CRISPRi-DR_results_nd_bsaf.txt"
  
  echo "Running CRISPRi-DR model for $drug with ND_bsAf base"
  transit cgi run_model "$fractional_abundances_nd_bsaf" "$crispri_results_nd_bsaf" -use_negatives
done

# Run add_rvnumber.py script for each CRISPRi-DR results file
for drug in $unique_drugs; do
  # ND_bsAf base
  input_file_nd_bsaf="ND_bsAf_results/${output_prefix}_${drug}_CRISPRi-DR_results_nd_bsaf.txt"
  output_file_nd_bsaf="ND_bsAf_results/${output_prefix}_${drug}_CRISPRi-DR_results_nd_bsaf_edited.txt"
  if [[ -f "$input_file_nd_bsaf" ]]; then
    python add_rvnumber.py -d "$input_file_nd_bsaf" -o "$output_file_nd_bsaf"
    echo "Processed $input_file_nd_bsaf -> $output_file_nd_bsaf"
  else
    echo "File $input_file_nd_bsaf not found"
  fi
done

echo "All tasks completed."
