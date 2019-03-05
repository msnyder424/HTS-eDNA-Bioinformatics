# HTS-eDNA-Bioinformatics
Scripts for the analysis of metabarcoding eDNA assay data from Illumina high-throughput sequencers

My PhD included the development and testing of environmental DNA assays for the detection of invasive species. These scripts are part of the bioinformatic pipeline to analyze those data.

MetaTrim.py trims primers and spacer inserts from Illumina paired end HTS reads. There are several methods for this program, which also works as a module that you can import and use various functions from. See the MetaTrimREADME

SeqTabToFasta.pl creates individual fasta files with amplicon sequence variants (ASVs) generated from the HTS reads merging program DADA2. The input (SeqTabNoChim.txt) is the sequence table output from DADA2 in tab delimited format.
