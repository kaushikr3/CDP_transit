#!/bin/bash

# Usage: ./run_transit.sh -i input_filelist.txt -r sgRNA_info.txt -m sample_metadata.txt -o output_prefix [-b base_name] -f

# Load necessary modules or environment
spack load /kinvl2y # r@4.2.2%gcc@8.2.0

# Parse arguments
bsaf=false
while getopts "i:r:m:o:b:f" opt; do
  case $opt in
    i) input_filelist="$OPTARG"
    ;;
    r) sgRNA_info="$OPTARG"
    ;;
    m) sample_metadata="$OPTARG"
    ;;
    o) output_prefix="$OPTARG"
    ;;
    b) base_name="$OPTARG"
    ;;
    f) bsaf=true
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
        exit 1
    ;;
  esac
done

# Check if input filelist is provided
if [ -z "$input_filelist" ];then
  echo "Input filelist is required. Usage: ./run_transit.sh -i input_filelist.txt -r sgRNA_info.txt -m sample_metadata.txt -o output_prefix [-b base_name] -f"
  exit 1
fi

# Check if sgRNA info file is provided
if [ -z "$sgRNA_info" ]; then
  echo "sgRNA info file is required. Usage: ./run_transit.sh -i input_filelist.txt -r sgRNA_info.txt -m sample_metadata.txt -o output_prefix [-b base_name] -f"
  exit 1
fi

# Check if sample metadata is provided
if [ -z "$sample_metadata" ]; then
  echo "Sample metadata is required. Usage: ./run_transit.sh -i input_filelist.txt -r sgRNA_info.txt -m sample_metadata.txt -o output_prefix [-b base_name] -f"
  exit 1
fi

# Check if output prefix is provided
if [ -z "$output_prefix" ]; then
  echo "Output prefix is required. Usage: ./run_transit.sh -i input_filelist.txt -r sgRNA_info.txt -m sample_metadata.txt -o output_prefix [-b base_name] -f"
  exit 1
fi

# Check if base name is provided when bsaf is true
if [ "$bsaf" = true ] && [ -z "$base_name" ]; then
  echo "Base name is required when -f is true. Usage: ./run_transit.sh -i input_filelist.txt -r sgRNA_info.txt -m sample_metadata.txt -o output_prefix -b base_name -f"
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
unique_drugs=$(python ~/CDP/scripts/extract_drugs.py "$sample_metadata")

if [ "$bsaf" = true ]; then
  # Create combined counts files for the specified base only
  base_files=($(printf "%s\n" "${count_files[@]}" | grep "$base_name" | sort))

  base_headers=()
  for file in "${base_files[@]}"; do
    base_headers+=($(basename "$file" .counts))
  done

  # Create combined counts file for the specified base
  if [ ${#base_files[@]} -gt 0 ]; then
    combined_counts_base="${output_prefix}_${base_name}_combined_counts.txt"
    echo "Creating combined counts for $base_name base"
    header_str_base=$(IFS=,; echo "${base_headers[*]}")
    count_files_str_base="${base_files[*]}"
    transit cgi create_combined_counts "$header_str_base" $count_files_str_base "$combined_counts_base"
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

      # Create combined headers and files for the specified base
      headers_base=("${headers[@]}" "${base_headers[@]}")
      combined_files_base=("${drug_files[@]}" "${base_files[@]}")
      combined_counts_base="${output_prefix}_${drug}_combined_counts_${base_name}.txt"
      if [ ${#base_files[@]} -gt 0 ]; then
        echo "Creating combined counts for $drug with $base_name base"
        header_str_base=$(IFS=,; echo "${headers_base[*]}")
        count_files_str_base="${combined_files_base[*]}"
        transit cgi create_combined_counts "$header_str_base" $count_files_str_base "$combined_counts_base"
      fi
    fi
  }

  # Generate create_combined_counts command for each drug
  for drug in $unique_drugs; do
    create_combined_counts_for_drug "$drug"
  done

  # Extract fractional abundances for each combined counts file
  for drug in $unique_drugs; do
    combined_counts_base="${output_prefix}_${drug}_combined_counts_${base_name}.txt"

    echo "Extracting fractional abundances for $drug with $base_name base"
    transit cgi extract_abund "$combined_counts_base" "$sample_metadata" "$base_name" "$sgRNA_info" Pred_logFC ORF "$drug" 1 "${output_prefix}_${drug}_frac_abundances_${base_name}.txt" -no_uninduced
  done

  # Run CRISPRi-DR model for each fractional abundance file
  for drug in $unique_drugs; do
    fractional_abundances_base="${output_prefix}_${drug}_frac_abundances_${base_name}.txt"

    crispri_results_base="${output_prefix}_${drug}_CRISPRi-DR_results_${base_name}.txt"

    echo "Running CRISPRi-DR model for $drug with $base_name base"
    transit cgi run_model "$fractional_abundances_base" "$crispri_results_base" -use_negatives
  done

  # Run add_rvnumber.py script for each CRISPRi-DR results file
  for drug in $unique_drugs; do
    # Specified base
    input_file_base="${base_name}_results/${output_prefix}_${drug}_CRISPRi-DR_results_${base_name}.txt"
    output_file_base="${base_name}_results/${output_prefix}_${drug}_CRISPRi-DR_results_${base_name}_edited.txt"
    if [[ -f "$input_file_base" ]]; then
      python ~/CDP/scripts/add_rvnumber.py -d "$input_file_base" -o "$output_file_base"
      echo "Processed $input_file_base -> $output_file_base"
    else
      echo "File $input_file_base not found"
    fi
  done
else
  # Create combined counts files for ND base only
  nd_files=($(printf "%s\n" "${count_files[@]}" | grep "ND_" | grep -v "ND_bsAf" | sort))

  nd_headers=()
  for file in "${nd_files[@]}"; do
    nd_headers+=($(basename "$file" .counts))
  done

  # Create combined counts file for ND
  if [ ${#nd_files[@]} -gt 0 ]; then
    combined_counts_nd="${output_prefix}_ND_combined_counts.txt"
    echo "Creating combined counts for ND base"
    header_str_nd=$(IFS=,; echo "${nd_headers[*]}")
    count_files_str_nd="${nd_files[*]}"
    transit cgi create_combined_counts "$header_str_nd" $count_files_str_nd "$combined_counts_nd"
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

      # Create combined headers and files for ND base
      headers_nd=("${headers[@]}" "${nd_headers[@]}")
      combined_files_nd=("${drug_files[@]}" "${nd_files[@]}")
      combined_counts_nd="${output_prefix}_${drug}_combined_counts_nd.txt"
      if [ ${#nd_files[@]} -gt 0 ]; then
        echo "Creating combined counts for $drug with ND base"
        header_str_nd=$(IFS=,; echo "${headers_nd[*]}")
        count_files_str_nd="${combined_files_nd[*]}"
        transit cgi create_combined_counts "$header_str_nd" $count_files_str_nd "$combined_counts_nd"
      fi
    fi
  }

  # Generate create_combined_counts command for each drug
  for drug in $unique_drugs; do
    create_combined_counts_for_drug "$drug"
  done

  # Extract fractional abundances for each combined counts file
  for drug in $unique_drugs; do
    combined_counts_nd="${output_prefix}_${drug}_combined_counts_nd.txt"

    echo "Extracting fractional abundances for $drug with ND base"
    transit cgi extract_abund "$combined_counts_nd" "$sample_metadata" ND "$sgRNA_info" Pred_logFC ORF "$drug" 1 "${output_prefix}_${drug}_frac_abundances_nd.txt" -no_uninduced
  done

  # Run CRISPRi-DR model for each fractional abundance file
  for drug in $unique_drugs; do
    fractional_abundances_nd="${output_prefix}_${drug}_frac_abundances_nd.txt"

    crispri_results_nd="${output_prefix}_${drug}_CRISPRi-DR_results_nd.txt"

    echo "Running CRISPRi-DR model for $drug with ND base"
    transit cgi run_model "$fractional_abundances_nd" "$crispri_results_nd" -use_negatives
  done

  # Run add_rvnumber.py script for each CRISPRi-DR results file
  for drug in $unique_drugs; do
    # ND base
    input_file_nd="ND_results/${output_prefix}_${drug}_CRISPRi-DR_results_nd.txt"
    output_file_nd="ND_results/${output_prefix}_${drug}_CRISPRi-DR_results_nd_edited.txt"
    if [[ -f "$input_file_nd" ]]; then
      python ~/CDP/scripts/add_rvnumber.py -d "$input_file_nd" -o "$output_file_nd"
      echo "Processed $input_file_nd -> $output_file_nd"
    else
      echo "File $input_file_nd not found"
    fi
  done
fi

echo "All tasks completed."
