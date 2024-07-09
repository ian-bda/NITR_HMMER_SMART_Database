# NITR_HMMER_SMART_Database

This is a universal version of Dustin Wcisel's perl scripts which can be found here: https://github.com/djwcisel/nilts, adapted to be able to handle any scaffold for the six frame translation and be able to search for NITR (Novel Immune Type Receptor) IG domains in Bichir.

## Step 1: Download SMART Database and HMMER update hmm files to be compatible
The SMART database was downloaded from [https://software.embl-em.de/](https://software.embl-em.de/software/18)

HMMER can be download from http://hmmer.org/

SMART.hmm file was made using hmmbuild.pl to update hmm files from SMART database to be compatible with HMMER by generating new .hmm profiles and merging into one file:

```
perl hmmbuild.pl
cat *.hmm > SMART.hmm
```

## Step 2: Grab Scaffolds:
Polypyerus bichir lapradei scaffolds were downloaded from: 

https://www.ncbi.nlm.nih.gov/nucleotide/JANRMT010028174.1?report%3Dgenbank%26log$%3Dnuclalign%26blast_rank%3D1%26RID%3D05P74N1E013&sa=D&source=docs&ust=1712765973402913&usg=AOvVaw301s9bpxW-7heRysVHYRGh

https://www.ncbi.nlm.nih.gov/nucleotide/JANRMT010062854.1?report=genbank&log$=nuclalign&blast_rank=1&RID=05P3DEH5013

https://www.ncbi.nlm.nih.gov/nucleotide/JANRMT010062807.1?report=genbank&log$=nuclalign&blast_rank=1&RID=05P74N1E013

## Step 3: Six Frame Translation (now universal!)

Create translations for DNA scaffolds

```
perl six_frame_univeral.pl JANRMT010028174.1.fa >> all_six_frames.fa
perl six_frame_univeral.pl JANRMT010062854.1.fa >> all_six_frames.fa 
perl six_frame_univeral.pl JANRMT010062807.1.fa >> all_six_frames.fa 
```

you want just one file containing all of the six frame transaltions for all of your scaffolds

## Step 5: Run hmmerscan on all_six_frames.fa file

```
hmmscan SMART.hmm all_six_frames.fa > all_six_frames_smart.out
```

## Step 6: Process HMMER Results for IG Domains

smart_hmmer_parser.pl searches the all_six_frames.out file for any ORFs that were significant for containing any variant of a SMART Immunoglobulin domain

```
perl smart_hmmer_parser.pl --scaffold scaffold.fasta > smart_parsed_hmmer_out.csv
```

replace "scaffold.fasta" with actual file. In this case: --scaffold JANRMT010028174.fasta, JANRMT010062854.fasta, and JANRMT010062807.fasta

# Find NITR domains

Now that I had a csv file of all the IG domains found in my scaffolds, I wanted to know which ones were NITRs. To do this I:

1. Grab protein sequences using:
   ```
   python Parse_HMMER_out.py smart_parsed_hmmer_out.csv
   ```
3. Create blast database from zebrafish protein assembly fasta file found here: https://www.ncbi.nlm.nih.gov/datasets/taxonomy/7955/ and blastp against ig_to_blast.fa --> blastp_out.txt
   
   ```
   makeblastdb -in zebrafish_protein.faa -dbtype prot -out zebrafish_protein_blast_db
   blastp -query ig_to_blast.fa -db zebrafish_protein_blast_db -out blastp_out.txt -max_target_seqs 5
   ```
6. Turn blastp_out.txt into csv file and merge with smart_parsed_hmmer_out.csv using nitr_blastp_parser.py 

    ```
   python nitr_blastp_parser.py blastp.txt smart_parsed_hmmer_out.csv
   ```

The final smart_parsed_hmmer_out.csv file should have the following columns: strand, scaffold, scaf_start, scaf_end, chr_start, chr_end, Ig-AA, Ig-nucl, BLASTP_hits which you can sort through and inspect manually for NITR domains.
