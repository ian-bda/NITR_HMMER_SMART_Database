#!/usr/bin/env python
# coding: utf-8

# ### Parse blastp out file into pandas dataframe

# In[42]:


import pandas as pd

# Initialize lists to store data
queries = []
alignments = []

# Read blastp output file
blastp_file = "path/to/zebrafish_protein_blastp_out"
with open(blastp_file) as f:
    lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("Query="):
            query = line.strip()
            queries.append(query)
            # Find the start index of the alignment section
            alignment_start_index = i + 6  # Skip lines until the start of alignment section
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

print(blastp_df)


# In[43]:


# Remove 'Query= ' from Query column 
blastp_df['Query'] = blastp_df['Query'].str.replace('Query= ', '', regex=True)
blastp_df


# In[44]:


# Remove 'Query= ' from Query column 
blastp_df['Query'] = blastp_df['Query'].str.replace('Query= ', '', regex=True)


# In[45]:


# Create scaf_start and scaf_end columns 
# Splitting the Query column
blastp_df[['Query', 'scaf_info']] = blastp_df['Query'].str.split('_', expand=True)
blastp_df[['scaf_start', 'scaf_end']] = blastp_df['scaf_info'].str.split('-', expand=True)

# Dropping not needed columns
blastp_df.drop(columns=['scaf_info'], inplace=True)
blastp_df = blastp_df.rename(columns={'Query': 'scaffold'})


# In[47]:


blastp_df


# In[59]:


# turn into csv file:
blastp_df.to_csv('blastp_out.csv',index = False)


# In[ ]:




