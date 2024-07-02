import pandas as pd
import sys

def extract_unique_drugs(sample_metadata_file):
    df = pd.read_csv(sample_metadata_file, sep='\t')
    df = df.iloc[:-6]
    unique_drugs = df['drug'].unique()
    return [str(drug) for drug in unique_drugs] 

if __name__ == "__main__":
    sample_metadata_file = sys.argv[1]
    unique_drugs = extract_unique_drugs(sample_metadata_file)
    for drug in unique_drugs:
        print(drug)

