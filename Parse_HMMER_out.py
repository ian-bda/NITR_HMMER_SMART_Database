#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd

def parse_hmmer_output(csv_file):
    try:
        # Read CSV into a DataFrame
        df = pd.read_csv(csv_file, sep='\t')
        
        # Check if the expected columns are present
        expected_columns = ['strand', 'scaffold', 'scaf_start', 'scaf_end', 'chr_start', 'chr_end', 'Ig-AA', 'Ig-nucl']
        for col in expected_columns:
            if col not in df.columns:
                print(f"Error: Column '{col}' not found in the CSV file.")
                return
        
        # Initialize an empty list to store sequences
        sequences = []
        
        # Iterate over rows in the DataFrame
        for index, row in df.iterrows():
            # Extract necessary columns
            scaffold = row['scaffold']
            scaf_start = row['scaf_start']
            scaf_end = row['scaf_end']
            sequence = row['Ig-AA']
            
            # Construct FASTA header
            header = f">{scaffold}_{scaf_start}-{scaf_end}\n"
            
            # Append header and sequence to the list
            sequences.append(header)
            sequences.append(sequence + "\n")
        
        # Write to FASTA file
        fasta_file = "ig_to_blast.fa"
        with open(fasta_file, "w") as f:
            f.writelines(sequences)
        
        print(f"The resulting {fasta_file} file has been created.")
    
    except FileNotFoundError:
        print(f"Error: File '{csv_file}' not found.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python Parse_HMMER_out.py <smart_parsed_hmmer_out.csv>")
    else:
        csv_file = sys.argv[1]
        parse_hmmer_output(csv_file)
