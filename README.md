# NITR_HMMER_SMART_Database

This is a universal version of Dustin Wcisel's perl scripts which can be found here: https://github.com/djwcisel/nilts, adapted to be able to handle any scaffold for the six frame translation and be able to search for NITR (Novel Immune Type Receptor) IG domains in Bichir.

## Step 1: Download SMART Database
The SMART database can be downloaded from https://software.embl-em.de/
````
perl hmmbuild.pl
cat *.hmm > SMART.hmm
hmmscan --cpu 12 SMART.hmm all_six_frames.fa > all_six_frames_smart.out
`````
## Step 2: Grab Scaffolds:
Polypyerus bichir lapradei NITR scaffolds were downloaded from: 

https://www.google.com/url?q=https://www.ncbi.nlm.nih.gov/nucleotide/JANRMT010028174.1?report%3Dgenbank%26log$%3Dnuclalign%26blast_rank%3D1%26RID%3D05P74N1E013&sa=D&source=docs&ust=1712765973402913&usg=AOvVaw301s9bpxW-7heRysVHYRGh

https://www.ncbi.nlm.nih.gov/nucleotide/JANRMT010062854.1?report=genbank&log$=nuclalign&blast_rank=1&RID=05P3DEH5013

https://www.ncbi.nlm.nih.gov/nucleotide/JANRMT010062807.1?report=genbank&log$=nuclalign&blast_rank=1&RID=05P74N1E013

## Step 3: Download and update HMMER to be compatible with SMART database
The SMART database was downloaded from [https://software.embl-em.de/](https://software.embl-em.de/software/18)

HMMER can be download from http://hmmer.org/

Use custom perl script to update hmm files from SMART database to be compatible with HMMER by generating new .hmm profiles and merging into one file:

```
perl hmmbuild.pl
cat *.hmm > SMART.hmm
hmmscan --cpu 12 SMART.hmm all_six_frames.fa > all_six_frames_smart.out
```


## Step 4: Six Frame Translation (now universal!)

Create translations for DNA scaffolds

```
perl six_frame_univeral.pl JANRMT010028174.1.fa >> all_six_frames.fa
perl six_frame_univeral.pl JANRMT010062854.1.fa >> all_six_frames.fa 
perl six_frame_univeral.pl JANRMT010062807.1.fa >> all_six_frames.fa 
```

# Step 5: 

smart_hmmer_parser.pl searches the all_six_frames.out file for any ORFs that were significant for containing any variant of a SMART Immunoglobulin domain

```
perl smart_hmmer_parser.pl > smart_parsed_hmmer_out.tsv
```

# Find NITR domains 

