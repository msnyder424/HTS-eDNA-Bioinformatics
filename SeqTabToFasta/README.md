# SeqTabToFasta.pl #

The R package Dada2 is a powerfull tool for removing sequencing error from high-throughput sequencing results. The output of this program is sequence table. The first row of this table is every sequence present in the samples processed, and every subsequent row is a histogram of the number of reads that were that sequence in each sample. This format is not useful for BLAST. SeqTabToFasta.pl takes a sequence table from Dada2 and creates individual fasta files for every sample in the input. Sequences in the fasta files will be named 1 through the total number of sequences and the number of reads that sequence represents separated by "|". The input file must be named "SeqTab.txt".

*Usage:*
`perl SeqTabToFasta.pl`
