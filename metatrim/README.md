# metatrim.py #

**by M. R. Snyder 2018.**

***Written for Python 3***

***metatrim*** is both a module with functions that can be called from a python interpreter and a stand alone script. This module will be packaged and available on pip upon completion of my PhD degree.

## Using *metatrim* ##

***metatrim*** requires one package that is not in the Python standard library. If this package is not installed in your environment it will be the first time you run ***metatrim***. Therefore, you must have the appropriate permissions to install Python packages the first time you use ***metatrim***.

***metatrim*** trims paired end metabarcoding reads returned from Illumina sequencing for subsequent merging in other programs such as Dada2, Unoise, or OBITools. ***metatrim*** works on the parent directory level. It should be run from the parent directory in which all of your subdirectories with Illumina sequencing results are stored. If you want to use individual functions from ***metatrim*** you can import it as a module. This allows you to trim individual sets of paired sequencing results files. ***metatrim*** handles zipped fastq.gz files. It will deposit your trimmed meta-barcoding results as unzipped fastq files into a new subdirectory named with your parent directory name and 'TrimmedFastqs'. If you have used primers with identical spacer inserts to those published in Klymus *et al.* 2017, *Plos One*, **12**(5): e0177643, you can use those inserts to remove instances of index hops. The most effective method of library prep for removing index hops is to combine every forward spacer with every reverse spacer. ***All inputs are case INsensitive! fastq.gz files must use the demultiplexed illumina naming convention (forward read file ends in R1.001.fastq.gz, reverse read file ends in R2_001.fastq.gz)***


#### *To trim primers from multiple Illumina HTS result files use **metatrim** as a stand alone script from the parent directory containing subdirectories for each sample with fastq.gz files:* ####
`python metatrim.py Primer_Set F_Seq R_Seq F_Err R_Err Length Spacers`  

Primer_Set: primer set name or "other" if it is not in the common primer set list. If a primer set name is entered, the next two argument variables will be ignored, but a value must be entered.  
F_Seq: Last N bases of forward primer. ≥8 is recomended.  
R_Seq: Last N bases of reverse primer. ≥8 is recomended.  
F_Err: N errors allowed in forward primer. ≥2 is recomended.  
R_Err: N errors allowed in reverse primer. ≥2 is recomended.   
Length: Length of marker OR primer set name to use a predefined length OR if length is variable, enter 0 and ***metatrim*** will search for the opposite primer in each forward and reverse read. If the primer is not found it will take the remainder of the sequence after the first primer is found. **WARNING:** if you choose to have ***metatrim*** search for the opposite primer in each read, you should ensure it is never found at any other location than is intended.  
Spacers: Are you using spacer inserts as published in Klymus et al. 2017, Plos One?: Yes(Y) or No(N). This option allows the removal of index hops!  


#### *To trim primers from a single sample use ***metatrim*** from a python interpreter:* ####
`MetaTrim(In_Forward, In_Reverse, Primer_Set, F_Seq R_Seq, F_Err, R_Err, Length Spacers)`

In_Forward: Path to forward ([...]_R1_001.fastq.gz) file.    
In_Reverse: Path to reverse ([...]_R2_001.fastq.gz) file.   
Primer_Set: Primer set name or "other" if it is not in the common primer set list. If a primer set name is entered, the next two variables will be ignored, but a value must be entered.  
F_Seq: Last N bases of forward primer. ≥8 is recomended.   
R_Seq: Last N bases of reverse primer. ≥8 is recomended.  
F_Err: N errors allowed in forward primer. ≥2 is recomended.  
R_Err: N errors allowed in reverse primer. ≥2 is recomended.    
Length: Length of marker OR primer set name to use a predefined length OR if length is variable, enter 0 and ***metatrim*** will search for the opposite primer in each forward and reverse read. If the primer is not found it will take the remainder of the sequence after the first primer is found. **WARNING:** if you choose to have ***metatrim*** search for the opposite primer in each read, you should ensure it is never found at any other location than is intended.  
Spacers: Are you using spacer inserts as published in Klymus et al. 2017, Plos One?: Yes(Y) or No(N). This option allows the removal of index hops!  

#### *To count spacers in a single fastq file* ####
*This feature is a function that cannot be called directly from \_\_main\_\_*  
`SpacerCount(In_fastq.gz)`  

In_fastq.gz: Path to fastq.gz file for which you want to count spacers.

#### *To view primer sets in the primer set list input:* ####
`python metatrim.py PrimerSets`  
or  
`PrintPrimerSets()`  


#### *To add a new primer set to the primer set list:* ####
`python metatrim.py New Primer_Set_Name F_Seq R_Seq Length`  
or  
`AddPrimer(Primer_Set_Name, F_Seq, R_Seq, Length)`  

Primer_Set_Name: Name your new primer set  
F_Seq: Last N bases of forward primer. ≥8 is recomended.   
R_Seq: Last N bases of reverse primer. ≥8 is recomended.  
Length: Length of marker OR primer set name to use a predefined length OR if length is variable, enter 0 and ***metatrim*** will search for the opposite primer in each forward and reverse read. If the primer is not found it will take the remainder of the sequence after the first primer is found. **WARNING:** if you choose to have ***metatrim*** search for the opposite primer in each read, you should ensure it is never found at any other location than is intended.  
*It is recommended that you input 8-10 bases of each primer. The sequence should end with the last 3' base of the primer.*  


#### *To remove a primer set from the list:* ####
`python metatrim.py Remove Primer_Set_Name`  
or  
`RemovePrimer(Primer_Set_Name)`  

Primer_Set_Name: Name of primer set to be removed

## metatrimExampleData ##
This folder contains some example data for users to familiarize themselves with the use of ***metatrim***. The two sub directories contain truncated gzipped fastq files. These samples can be trimmed with the *CYPPART* primer set with or without removing incorrect spacer inserts.  
### *Example usage* ###  
`python metatrim.py Primer_Set F_Seq R_Seq F_Err R_Err Length Spacers`  

Primer_Set: 'CyPpArT'. Note: ***metatrim*** is case INsensitive.  
F_Seq: Enter any value. ***metatrim*** will use the CYPPART primer set.  
R_Seq: Enter any value. ***metatrim*** will use the CYPPART primer set.  
F_Err: N errors allowed in forward primer. ≥2 is recomended.  
R_Err: N errors allowed in reverse primer. ≥2 is recomended.   
Length: ***metatrim*** will use the CYPPART primer set length unless you enter a value. You can enter any length, but a value of 0 will result in ***metatrim*** searching for both primers in every sequence.  
Spacers: Enter 'Y' or 'N'.    
