#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd

def parse_blastp(blastp_file):
    # Initialize lists to store data
    queries = []
    alignments = []

    # Read the BLASTP output file
    with open(blastp_file) as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            line = lines[i]
            if line.startswith("Query="):
                query = line.strip()
                queries.append(query)
                # Find the start index of the alignment section
                alignment_start_index = i + 6  # Skip lines until the start of the alignment section
                # Find the end index of the alignment section
                alignment_end_index = alignment_start_index
                for j in range(alignment_start_index, len(lines)):
                    if lines[j].strip() == "":
                        alignment_end_index = j
                        break
                # Extract alignment section
                alignment = "".join(lines[alignment_start_index:alignment_end_index])
                alignments.append(alignment)
                # Move to the next query
                i = alignment_end_index + 4
            else:
                i += 1

    # Create DataFrame
    blastp_df = pd.DataFrame({
        "Query": queries,
        "Alignment": alignments
    })

    # Remove 'Query= ' from Query column
    blastp_df['Query'] = blastp_df['Query'].str.replace('Query= ', '', regex=True)

    # Split the Query column to extract scaffold information
    blastp_df[['scaffold', 'scaf_start', 'scaf_end']] = blastp_df['Query'].str.extract(r'(.+?)_(\d+)-(\d+)')
    
    return blastp_df

def main(blastp_file, hmmer_file):
    blastp_df = parse_blastp(blastp_file)
    
    # Load the HMMER parsed output with tab delimiter
    hmmer_df = pd.read_csv(hmmer_file, delimiter='\t')

    # Print column names for debugging
    print("HMMER DataFrame columns:", hmmer_df.columns)
    
    # Initialize a new column for BLASTP hits
    hmmer_df['BLASTP_hits'] = ""

    # Iterate through HMMER DataFrame to match and update BLASTP results
    for idx, row in hmmer_df.iterrows():
        scaffold = row['scaffold']
        scaf_start = str(row['scaf_start'])
        scaf_end = str(row['scaf_end'])
        
        # Find matching BLASTP entries
        matches = blastp_df[(blastp_df['scaffold'] == scaffold) & 
                            (blastp_df['scaf_start'] == scaf_start) & 
                            (blastp_df['scaf_end'] == scaf_end)]
        
        if not matches.empty:
            # Combine the top 5 alignments
            top_hits = matches['Alignment'].head(5).tolist()
            hmmer_df.at[idx, 'BLASTP_hits'] = "\n".join(top_hits)

    # Save the updated DataFrame back to CSV
    hmmer_df.to_csv(hmmer_file, index=False, sep='\t')

    print(f"CSV file '{hmmer_file}' has been updated successfully.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python parse_blastp.py <blastp_output_file> <hmmer_output_file>")
    else:
        blastp_file = sys.argv[1]
        hmmer_file = sys.argv[2]
        main(blastp_file, hmmer_file)
