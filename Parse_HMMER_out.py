#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from Bio import SeqIO


# In[2]:


df = pd.read_csv('path/to/smart_parsed_hmmer_out.csv')


# In[3]:


df


# In[4]:


sequences = df['Ig-AA']
scaf_start = df['scaf_start']
scaf_end = df['scaf_end']
scaffold = df['scaffold']


# In[5]:


fasta_file = "ig_to_blast.fa"
with open(fasta_file, "w") as f:
    for index, sequence in enumerate(sequences):
        # Construct FASTA header
        header = f">{scaffold[index]}_{scaf_start[index]}-{scaf_end[index]}\n"
        # Write header and sequence to file
        f.write(header)
        f.write(sequence + "\n")


# In[ ]:




