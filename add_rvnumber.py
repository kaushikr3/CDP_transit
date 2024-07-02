import pandas as pd
import argparse

def process_dataframe(input_file, output_file):
    # Read the file to find the last comment line
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Find the last comment line
    for i in range(len(lines) - 1, -1, -1):
        if lines[i].startswith('#'):
            header_line = lines[i].strip('#').strip()
            break

    # Create a new temporary file with the proper header
    temp_file_path = 'temp_file.txt'
    with open(temp_file_path, 'w') as temp_file:
        temp_file.write(header_line + '\n')
        for line in lines[i+1:]:
            temp_file.write(line)

    # Load the dataset
    df = pd.read_csv(temp_file_path, sep='\t')

    # Display the first few rows of the DataFrame
    print("Initial DataFrame:")
    print(df.head())

    # Add Rv_number column and rearrange the df
    df['Rv_number'] = df.Orf.str.split(':').str[0].replace("RVBD", "rv", regex=True)
    df.insert(3, "Rv_number", df.pop("Rv_number"))

    # Save the DataFrame to a tab-separated file
    df.to_csv(output_file, sep='\t', index=False)
    print(f"\nDataFrame saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Process a DataFrame to add a new column 'Rv_number'.")
    parser.add_argument('-d', '--data', required=True, help='Path to the input data file')
    parser.add_argument('-o', '--output', required=True, help='Path to save the output data file')

    args = parser.parse_args()
    process_dataframe(args.data, args.output)

if __name__ == "__main__":
    main()
