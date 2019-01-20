# HTS-eDNA-Bioinformatics
Scripts for the analysis of metabarcoding eDNA assay data from Illumina high-throughput sequencers

My PhD included the development and testing of environmental DNA assays for the detection of invasive species. These scripts are part of the bioinformatic pipeline to analyze those data.

MetaTrim.py trims primers and spacer inserts from Illumina paired end HTS reads. There are several methods for this program. See the MetaTrimREADME

SeqTabToFasta.pl creates individual fasta files with amplicon sequence variants (ASVs) generated from the HTS reads merging program DADA2. The input (SeqTabNoChim.txt) is the sequence table output from DADA2 in tab delimited format.

SummarizeBLAST.pl summarizes BLAST resulst. It combines all hits with the lowest e value for a single query. Hits for the same species or groups of species are counted and the output has the proportion of reads in a sample that were of that sequence. See the SummarizeBLASTREADME.

CombineMarkers.py compares blast results from multiple eDNA assays. It identifies hits that were below an error cutoff (default is 0.1% of sequences, but can use errors caculated from positive controls per sample as an input file or any user defined cuttoff including 0%). It will optionally remove these hits, or if a species was found in multiple assays it will retain those hits, regardless of the frequency of the sequence.
