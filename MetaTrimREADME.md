MetaTrim.py

by M. R. Snyder 2018. 

Written for Python 3

MetaTrim trims paired end metabarcoding reads returned from Illumina sequencing for subsequent merging in other programs such as Dada2, Unoise, or OBITools. MetaTrim works on the parent directory level. It should be run from the parent directory in which all of your subdirectories with Illumina sequencing results are stored. If you want to use individual functions from MetaTrim you can import it as a module. This allows you to trim individual sets of paired sequencing results files. MetaTrim handles zipped fastq.gz files. It will deposit your trimmed meta-barcoding results as unzipped fastq files into a new subdirectory named with your parent directory name and 'TrimmedFastqs'. If you have used primers with identical spacer inserts to those published in Klymus et al. 2017 in Plos One, you can use those inserts to remove instances of index hops. The most effective method of library prep for removing index hops is to combine every forward spacer with every reverse spacer.

To use MetaTrim as a stand alone script input the following argument variables: 
1. primer set name or "other" if it is not in the common primer set list. If a primer set name is entered, argument variables 2 & 3 will be ignored, but a value must be entered.
2. Last N bases of forward primer 
3. Last N bases of reverse primer
4. N errors allowed in forward primer
5. N errors allowed in reverse primer 
6. Length of marker OR primer set name to use a predefined length OR if length is variable, enter 0 and MetaTrim will search for the opposite primer in each forward and reverse read. If the primer is not found it will take the remainder of the sequence after the first primer is found. WARNING: if you choose to have MetaTrim search for the opposite primer in each read, you should ensure it is never found at any other location than is intended.
7. 'Are you using spacer inserts as published in Klymus et al. 2017, Plos One?: Y/N'

To add a new primer set to the primer set list input the following argument vaiables: 1. 'New' in primer set. 2. Primer_Set Name, 3. forward sequence, 4. reverse sequence, 5. Length of the marker or input 0 (zero) if you want MetaTrim to always search for the opposite primer in each read. 

To remove a primer set from the list input the following argument variables: 1. Remove. 2. Primer_Set_Name
