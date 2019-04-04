#!/usr/bin/env python
"""
MetaTrim

by M. R. Snyder 2018. 
Written for Python 3

MetaTrim trims paired end metabarcoding reads returned from Illumina sequencing for subsequent 
merging in other programs such as Dada2, Unoise, or OBITools. metatrim.py works on the parent 
directory level. It should be run from the parent directory in which all of your subdirectories 
with Illumina sequencing results are stored. If you want to use individual functions from MetaTrim 
you can import it as a module. This allows you to trim individual sets of paired sequencing results 
files. MetaTrim handles zipped fastq.gz files from a python interpreter. It will deposit your trimmed 
meta-barcoding results as unzipped fastq files into a new subdirectory named with your parent directory 
name and 'TrimmedFastqs'. If you have used primers with identical spacer inserts to those published 
in Klymus et al. 2017 in Plos One "Environmental DNA (eDNA) metabarcoding assays to detect invasive invertebrate 
speciesin the Great Lakes" (https://doi.org/10.1371/journal.pone.0177643), you can use those inserts to remove 
instances of index hops. The most effective method of library prep for removing index hops is to combine every 
forward spacer with every reverse spacer.

Refer to the MetaTrim README.md for usage:
https://github.com/msnyder424/HTS-eDNA-Bioinformatics/blob/master/metatrim/README.md

Example data is available at:
https://github.com/msnyder424/HTS-eDNA-Bioinformatics/tree/master/metatrim
"""

import re
from itertools import permutations
import sys
import gzip
import os
import operator
from datetime import datetime

if __name__ != '__main__':
    print("MetaTrim takes the following variables:\n\
        1. primer set name or 'other' if it is not in the \
        common primer set list. If a primer set name is entered, argument variables 2 & 3 will be ignored, \
        but a value must be entered.\n\
        2. Last N bases of forward primer (recomended >= 8)\n\
        3. Last N bases of reverse primer (recomended >= 8)\n\
        4. N errors allowed in forward primer (recomended >=1)\n\
        5. N errors allowed in reverse primer (recomended >=1)\n\
        6. Length of marker. If you wish to use the length associated with a primer set, input the primer \
        set name. If length is variable, enter 0 and MetaTrim will search for the opposite primer in each \
        forward and reverse read. If the primer is not found it will take the remainder of the sequence \
        after the first primer is found. \
        WARNING: if your marker is much shorter than the read length, you should ensure that the oposite primer \
        sequence is not found in the Illumina sequencing primer or in any region after it's intended \
        location.\n\
        7. Are you using spacer inserts as published in Klymus et al. 2018, Plos One?: Y/N\n\n\
        To add a new primer set to the primer set list use AddPrimer and input the following vaiables: \
        1. Primer_Set Name, \
        2. Forward sequence, \
        3. Reverse sequence, \
        4. Length of the marker, primer set name (to use the stored length), or input 0 if you want MetaTrim \
        to always search for the opposite primer in each read. If the primer is not found it will take the \
        remainder of the sequence after the first primer is found. \
        WARNING: if you choose to have MetaTrim search for the oposite primer in each read, you should \
        ensure it is never found at any other location than is intended. This is most problematic when \
        your target can be very short and the oposite primer is found within the Illumina sequencing primer.\n\
        To remove a primer set from the list use RemovePrimer and imput the following variables: \
        1. Primer_Set_Name")
"""
PrimerSets: this dictionary contains commonly used primer sets. Each key is the primer set name.
The value is a list with structure ['ForwardSequence', 'ReverseSequence', 'MarkerLength']. See MetaTrimREADME for more info.
Primers in the list that comes with MetaTrim.py are from Snyder et al. 2019 
"Invasive species in bait and pond stores: metabarcoding environmental DNA assays and angler, 
retailer, and manager implications"
"""
PrimerSets = {
	'ACTSPART-2': ['TNACNTTCCGN','CNCCAATTCAN','53'],
    'MIFISHPART-2': ['TCGTGCCAGCN','TCCCAGTTTGN','0'],
    'ACTLPART-2': ['GMTCHATYCCN','CNCCAATTCAN','150'],
    'GOBIPART-2': ['TWAAAATYGCN','ACRTCWCGRCN','165'],
    'CYPPART-2': ['CYCTHCTAGGN','CYCCRTTRGCN','134'],
    'CYPCOMPLETE': ['TGATGAAAYTTYGGMTCYCTHCTAGG','AARAAGAATGATGCYCCRTTRGC','136'],
    'GOBICOMPLETE': ['AACVCAYCCVCTVCTWAAAATYGC','AGYCANCCRAARTTWACRTCWCGRC','165'],
    'MIFISHPART': ['TCGTGCCAGC','TCCCAGTTTG','0'],
    'ACTLPART': ['GMTCHATYCC','CNCCAATTCA','152'],
    'GOBIPART': ['TWAAAATYGC','ACRTCWCGRC','167'],
    'CYPPART': ['CYCTHCTAGG','CYCCRTTRGC','136'],
    'CYPCOMPLETE': ['TGATGAAAYTTYGGMTCYCTHCTAGG','AARAAGAATGATGCYCCRTTRGC','136'],
    'GOBICOMPLETE': ['AACVCAYCCVCTVCTWAAAATYGC','AGYCANCCRAARTTWACRTCWCGRC','165']
    }

#Finction to print primer sets in the primer set list
def PrintPrimerSets():
    print ('Primer Set\tFor Seq\tRev Seq\tLength')
    for keys in PrimerSets:
        print ('%s\t%s\t%s\t%s' % (keys, PrimerSets[keys][0], PrimerSets[keys][1], PrimerSets[keys][2]))

#Add a new primer to the primer set list
def AddPrimer(PrimerName, PF, PR, LenMarker):
    if PrimerName in PrimerSets:
        print (PrimerName, "is already in the primer set list!")
        exit()
    else:
        NewOrNot = input("Would you like to rename this version of MetaTrim? Yes(Y), No(N), or Cancel (X)").upper()
        if NewOrNot == 'Y':
            NewPyFileName = input("Name the new MetaTrim program file. Don't forget .py!")
        elif NewOrNot == 'N':
            NewPyFileName = 'metatrim.py'
        elif NewOrNot == 'X':
            exit()
        print ('Adding primer set:', PrimerName, 'F seq:', PF, 'R seq:', PR, 'Target length:', \
            LenMarker, 'to python file', NewPyFileName)
        with open('metatrim.py', 'r') as pyfileOld:
            lines = pyfileOld.readlines()
        for indx, val in enumerate(lines):
            if re.search('PrimerSets \= \{', val):
                NewPLineIndx = indx + 1
                lines.insert(NewPLineIndx, "\t'"+PrimerName.upper()+"': ['"+PF+"','"+PR+"','"+str(LenMarker)+"'],\n")
                break
        with open(NewPyFileName, 'w') as pyfileNew:
            for i in lines:
                pyfileNew.write(i)

#Remove a primer from the primer set list
def RemovePrimer(PrimerName):
    if not PrimerName.upper() in PrimerSets:
        print (PrimerName.upper(), "is not in the primer set list!")
        exit()
    else:
        NewOrNot = input("Would you like to rename this version of MetaTrim?(Y/N)").upper()
        if NewOrNot == 'Y':
            NewPyFileName = input("Name the new MetaTrim program file. Don't forget .py!")
        else:
            NewPyFileName = 'metatrim.py'
        PrimerToRemove = "\t'"+PrimerName.upper()+"':"
        print ('Removing primer set', PrimerName.upper(), 'from python file', NewPyFileName)
        with open('metatrim.py', 'r') as pyfileOld:
            lines = pyfileOld.readlines()
            for i in lines:
                if re.match(PrimerToRemove, i):
                    lines.remove(i)
                if re.match('#EndPrimerSetList', i):
                    break
        with open(NewPyFileName, 'w') as pyfileNew:
            for i in lines:
                pyfileNew.write(i)

#Function for counting correct spacer insert
FSpacers = {'e': '(TCCT|[ATCGN]CCT|T[ATCGN]CT|TC[ATCGN]T|TCC[ATCGN])', 
            'f': '(ATGC|[ATCGN]TGC|A[ATCGN]GC|AT[ATCGN]C|ATG[ATCGN])', 
            'g': '(CGAG|[ATCGN]GAG|C[ATCGN]AG|CG[ATCGN]G|CGA[ATCGN])', 
            'h': '(GATA|[ATCGN]ATA|G[ATCGN]TA|GA[ATCGN]A|GAT[ATCGN])'
            }
RSpacers = {'e': '(CGTA|[ATCGN]GTA|C[ATCGN]TA|CG[ATCGN]A|CGT[ATCGN])', 
            'f': '(TCAC|[ATCGN]CAC|T[ATCGN]AC|TC[ATCGN]C|TCA[ATCGN])', 
            'g': '(GAGT|[ATCGN]AGT|G[ATCGN]GT|GA[ATCGN]T|GAG[ATCGN])', 
            'h': '(ATCG|[ATCGN]TCG|A[ATCGN]CG|AT[ATCGN]G|ATC[ATCGN])'
            }

def SpacerCount(InFastq):
    Count = {'e': 0, 'f': 0, 'g': 0, 'h': 0}
    readbuffer = []
    global reads
    global MaxS
    reads = 0
    if re.search('R1_001\.fastq\.gz', InFastq):
        ReadDirection = 'forward'
        Spacer = FSpacers
    elif re.search('R2_001\.fastq\.gz', InFastq):
        ReadDirection = 'reverse'
        Spacer = RSpacers
    print ('Removing index hops based on incorrect spacer insert in sample', \
        InFastq[0:re.search('_', InFastq).start()], ReadDirection)
    with gzip.open(InFastq) as infile:
        for line in infile:
            readbuffer.append(line.strip().decode('ascii'))
            if len(readbuffer) == 4:
                reads += 1
                for keys in Spacer:
                    if re.match(Spacer[keys], readbuffer[1]):
                        Count[keys] += 1
                readbuffer = []
            if reads == 1000:
                break
    print ('Spacer counts in first 1K sequences:')
    for keys in Count:
        print (keys, ':', Count[keys])
    MaxS = max(Count.items(), key=operator.itemgetter(1))[0]
    print (InFastq[0:re.search('_', InFastq).start()], ReadDirection, 'has spacer', MaxS)

#Reverse compliment a sequence
def RevComp(Seq):
    global SeqRC
    trans = str.maketrans('ATGCN', 'TACGN')
    SeqRC = Seq[::-1].translate(trans)


#Create degenerate primer regexes
IUPACAmb = {'R' : '[AG]', 
            'Y' : '[CT]', 
            'S' : '[GC]', 
            'W' : '[AT]', 
            'K' : '[GT]', 
            'M' : '[AC]', 
            'B' : '[CGT]', 
            'D' : '[AGT]', 
            'H' : '[ACT]', 
            'V' : '[ACG]', 
            'N' : '[ATCG]'
            }
def DegPrimers (PrimDict, ErrDict):
    global DegPrimerDict
    DegPrimerDict = {}
    print ('Making degenerate primer regex.')
    for keys in PrimDict:
        DegPrimers = []
        PrimCurr = list(PrimDict[keys])
        PrimDegDict = {}
        if ErrDict[keys] > 0:
            for err in range(0, ErrDict[keys]+1):
                PrimDeg = []
                for bp in range(len(PrimDict[keys])):
                    if bp < err:
                        PrimDeg.append(0)
                    else:
                        PrimDeg.append(1)
                Perms = permutations(PrimDeg)
                for i in Perms:
                    if i in PrimDegDict:
                        continue
                    else:
                        PrimDegDict[i] = 1
            for i in PrimDegDict:
                DegCurr = list(i)
                for n in range(len(DegCurr)):
                    if DegCurr[n] == 0:
                        DegCurr[n] = '[ACTGN]'
                    else:
                        if PrimCurr[n] in IUPACAmb:  
                            DegCurr[n] = IUPACAmb[PrimCurr[n]]
                        else:
                            DegCurr[n] = PrimCurr[n]
                DegPrimers.append(''.join(DegCurr))
                DegPrim = '|'.join(DegPrimers)
            DegPrimerDict[keys] = '('+DegPrim+')'
        else:
            DegCurr = []
            for n in PrimCurr:
                if n in IUPACAmb:  
                    DegCurr.append(IUPACAmb[n])
                else:
                    DegCurr.append(n)
            DegPrimers.append(''.join(DegCurr))
            DegPrim = '|'.join(DegPrimers)
            DegPrimerDict[keys] = '('+DegPrim+')'

#Function to trim primers
def TrimPrimers (Primer1,Primer2):
    global name
    global TargetStart
    global TargetEnd
    global TrimmedSeq
    global Seq
    global Qual
    if re.search(DegPrimerDict[Primer1], readbuffer[1]):
        name = readbuffer[0].split()
        TargetStart = re.search(DegPrimerDict[Primer1], readbuffer[1]).end()
        if Length == 0:
            RevComp(readbuffer[1])
            if re.search(DegPrimerDict[Primer2], SeqRC):
                TargetEnd = len(readbuffer[1]) - re.search(DegPrimerDict[Primer2], SeqRC).end()
            else:
                TargetEnd = len(readbuffer[1])
        else:
            TargetEnd = TargetStart + int(Length)
        Seq = readbuffer[1][TargetStart:TargetEnd]
        Qual = readbuffer[3][TargetStart:TargetEnd]
    else:
        Seq = 'X'
        name = readbuffer[0].split()

def MetaTrim(InForward, InReverse, PrimerSet, PF, PR, ErrF, ErrR, TargetLen, Spacers):
    start = datetime.now().time()
    PrimerSet = PrimerSet.upper()
    global Length
    global PrimerDict
    global ErrorDict
    global cutoff
    if type(TargetLen) == int:
        Length = int(TargetLen)   
    else:
        Length = TargetLen.upper()      

    #Define the primer set if not in the list
    if PrimerSet == 'OTHER':
        print ('Primer set is ', PrimerSet,'. F: ', PF.upper(), 'R: ', PR.upper())
        PrimerDict = {'pF': PF.upper(), 'pR': PR.upper()}
        ErrorDict = {'pF': int(ErrF), 'pR': int(ErrR)}
    elif PrimerSet in PrimerSets:
        print ('Primer set is ', PrimerSet)
        PrimerDict = {'pF': PrimerSets[PrimerSet][0], 'pR': PrimerSets[PrimerSet][1]}
        ErrorDict = {'pF': int(ErrF), 'pR': int(ErrR)}
        if Length in PrimerSets:
            Length = int(PrimerSets[PrimerSet][2])
    elif PrimerSet not in PrimerSets:
        print (PrimerSet.upper(),"is not in the Primer Sets List!")
        if __name__ == '__main__':
            exit()

    #Determine length of sequences to remove
    try:
        cutoff
    except NameError:
        if Length > 125:
            cutoff = 100
        elif Length != 0:
            cutoff = Length * 0.75
        elif Length == 0:
            cutoff = int(input('What size of sequences should be considered primer dimer and removed?'))

    #Create Ddegenerate primer regex if it does not already exist
    try:
        DegPrimerDict
    except NameError:
        DegPrimers (PrimerDict, ErrorDict)

    cwd = os.getcwd()
    basenm = os.path.basename(cwd)
    outsumname = basenm+'TrimSummary.txt'

    #Create summary files

    if __name__ != '__main__':
        outsum = open(outsumname, "w")
        if Spacers == 'Y':
            outsum.write('Sample\tReads\tF Seqs Trimmed\tR Seqs Trimmed\tSeqs F & R Trimmed\tShort Seqs\t\
                Seqs w/ Correct Spacer Combo\tSeqs w/o Correct Spacer Combo\n')
        else:
            outsum.write('Sample\tReads\tF Seqs Trimmed\tR Seqs Trimmed\tSeqs F & R Trimmed\tShort Seqs\n')
    else:
        outsum = open(outsumname, "a")

    #Create results directory
    ResDirName = basenm+'TrimmedFastqs'
    try:
        os.mkdir(ResDirName)
        print("Directory " , ResDirName ,  " Created ") 
    except FileExistsError:
        pass

    nFSeqs = 0
    FNames = {}
    FSeqs = {}
    FQuals = {}
    CorSpacer = {}
    RSeqs = 0
    FinalSeqs = 0
    Shorts = 0
    sample = InForward[0:re.search('_', InForward).start()]
    if re.search('R1_001\.fastq\.gz', InForward):
        reads = 0
        Foutfilenameandpath = ResDirName+'/'+sample+'_R1_001.fastq'
        Foutfile = open(Foutfilenameandpath, "w")
        global readbuffer
        readbuffer = []
        print ('Processing:', InForward)
        if Spacers.upper() == 'Y':
            #count spacer inserts in firs 1K sequences
            SpacerCount(InForward)
            FMax = MaxS
            print ('Trimming primers sample', sample, 'forward.')
            with gzip.open(InForward) as infile:
                for line in infile:
                    readbuffer.append(line.strip().decode('ascii'))
                    if len(readbuffer) == 4:
                        reads += 1
                        if reads % 10000 == 0:
                            print ('Read:', reads, end='\r')
                        if re.match(FSpacers[FMax], readbuffer[1]):
                            TrimPrimers('pF','pR')
                            if re.match('[ATCGN]', Seq):
                                nFSeqs += 1
                                if len(Seq) > cutoff:
                                    FNames[name[0]] = readbuffer[0]
                                    FSeqs[name[0]] = Seq
                                    FQuals[name[0]] = Qual
                                    CorSpacer[name[0]] = 1
                        readbuffer = []
            print(reads, 'total raw reads.')
        else:
            with gzip.open(InForward) as infile:
                print ('Trimming primers sample', sample, 'forward.')
                for line in infile:
                    LineCurr = line.strip().decode('ascii')
                    readbuffer.append(LineCurr)
                    if len(readbuffer) == 4:
                        reads += 1
                        if reads % 10000 == 0:
                            print ('Read:', reads, end='\r')
                        TrimPrimers('pF','pR')
                        if re.match('[ATCGN]', Seq):
                            nFSeqs += 1
                            if len(Seq) > cutoff:
                                FNames[name[0]] = readbuffer[0]
                                FSeqs[name[0]] = Seq
                                FQuals[name[0]] = Qual
                        readbuffer = []
                print(reads, 'total raw reads.')
    else:
        print ("MetaTrim only works on paired zipped fastq files!")
    if re.search('R2_001\.fastq\.gz', InReverse):
        Routfilenameandpath = ResDirName+'/'+sample+'_R2_001.fastq'
        Routfile = open(Routfilenameandpath, "w")
        readbuffer = []
        reads = 0
        print ('Processing:', InReverse)
        if Spacers.upper() == 'Y':
            #count spacer inserts in firs 1K sequences
            SpacerCount(InReverse)
            RMax = MaxS
            print ('Trimming primers sample', sample, 'forward.')
            with gzip.open(InReverse) as infile:
                for line in infile:
                    LineCurr = line.strip().decode('ascii')
                    readbuffer.append(LineCurr)
                    if len(readbuffer) == 4:
                        reads += 1
                        if reads % 10000 == 0:
                            print ('Read:', reads, end='\r')
                        if re.match(RSpacers[RMax], readbuffer[1]):
                            TrimPrimers('pR','pF')
                            if re.match('[ATCGN]', Seq):
                                RSeqs += 1
                                if name[0] in FSeqs:
                                    if len(Seq) > cutoff:
                                        Foutfile.write("%s\n%s\n+\n%s\n" % (FNames[name[0]],FSeqs[name[0]],FQuals[name[0]]))
                                        Routfile.write("%s\n%s\n+\n%s\n" % (readbuffer[0],Seq,Qual))
                                        FinalSeqs += 1
                                    else:
                                        Shorts += 1
                        else:
                            CorSpacer[name[0]] = 0
                        readbuffer = []
                print(reads, 'total raw reads.')
        else:
            with gzip.open(InReverse) as infile:
                print ('Trimming primers sample', sample, 'reverse.')
                for line in infile:
                    LineCurr = line.strip().decode('ascii')
                    readbuffer.append(LineCurr)
                    if len(readbuffer) == 4:
                        reads += 1
                        if reads % 10000 == 0:
                            print ('Read:', reads, end='\r')
                        TrimPrimers('pR','pF')
                        if re.match('[ATCGN]', Seq):
                            RSeqs += 1
                            if name[0] in FSeqs:
                                if len(Seq) > cutoff:
                                    Foutfile.write("%s\n%s\n+\n%s\n" % (FNames[name[0]],FSeqs[name[0]],FQuals[name[0]]))
                                    Routfile.write("%s\n%s\n+\n%s\n" % (readbuffer[0],Seq,Qual))
                                    FinalSeqs += 1
                                else:
                                    Shorts += 1
                        readbuffer = []
                print(reads, 'total raw reads.')
        Foutfile.close()
        Routfile.close()
        if Spacers == 'Y':
            CorrSpacerCount = sum(CorSpacer.values())
            IncorrSpacerCount = reads - CorrSpacerCount
            outsum.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sample, reads, nFSeqs, RSeqs, \
                FinalSeqs, Shorts, CorrSpacerCount, IncorrSpacerCount))
        else:
            outsum.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (sample, reads, nFSeqs, RSeqs, FinalSeqs, Shorts))
    else:
        print ("MetaTrim only works on paired zipped fastq files!")
    end = datetime.now().time()
    if __name__ != '__main__':
        print ('MetaTrim start:', start, '\n', 'MetaTrim end:', end) 
        outsum.close()


if __name__ == "__main__":
    try: 
        sys.argv[1]
    except IndexError:
        print("To use input the following argument variables:\n\
            1. primer set name or 'other' if it is not in the common primer set list. \
            If a primer set name is entered, argument variables 2 & 3 will \
            be ignored, but a value must be entered.\n\
            2. Last N bases of forward primer (recomended >= 8)\n\
            3. Last N bases of reverse primer (recomended >= 8)\n\
            4. N errors allowed in forward primer (recomended >=1)\n\
            5. N errors allowed in reverse primer (recomended >=1)\n6. Length of marker. If you wish to use \
            the length associated with a primer set, input the primer set name. If length is variable, \
            enter 0 and MetaTrim will search for the opposite primer in each forward and reverse read. \
            If the primer is not found it will take the remainder of the sequence after the first primer \
            is found. WARNING: if your marker is much shorter than the read length, you should ensure that the \
            oposite primer sequence is not found in the Illumina sequencing primer or in any region after it's \
            intended location.\n\
            7. Are you using spacer inserts as published in Klymus et al. 2018, Plos One?: Y/N\n\n\
            To add a new primer set to the primer set list input the following argument vaiables: \
            1. 'New' in primer set. \
            2. Primer_Set Name, \
            3. Forward sequence, \
            4. Reverse sequence, \
            5. Length of marker OR primer set name to use a predefined length OR if length is variable, enter 0 \
            and MetaTrim will search for the opposite primer in each forward and reverse read. If the primer is \
            not found it will take the remainder of the sequence after the first primer is found. \
            WARNING: if you choose to have MetaTrim search for the oposite primer in each read, you should \
            ensure it is never found at any other location than is intended. This is most problematic when \
            your target can be very short and the oposite primer is found within the Illumina sequencing primer. \
            It is recommended that you input >8 bases of each primer. the sequence should end with the last \
            3' base of the primer.\n\
            To remove a primer set from the list input the following argument variables: \
            1. Remove. \
            2. Primer_Set_Name")
    else:
        start = datetime.now().time()
        #Add a new primer to the primer set list
        if sys.argv[1].upper() == 'PRIMERSETS':
            PrintPrimerSets()
            exit()
        elif sys.argv[1].upper() == 'NEW':
            if len(sys.argv) == 6:
                AddPrimer(sys.argv[2].upper(), sys.argv[3].upper(), sys.argv[4].upper(), sys.argv[5])
                exit()
            else:
                print ("To add a new primer set to the primer set list input the following argument vaiables: \
                    1. 'New'. 2. Primer_Set Name, \
                    3. Forward sequence, 4. Reverse sequence, \
                    5. Length of marker OR primer set name to use a predefined length OR if length is variable, \
                    enter 0 and MetaTrim will search for the opposite primer in each forward and reverse read. \
                    If the primer is not found it will take the remainder of the sequence after the first primer \
                    is found. WARNING: if you choose to have MetaTrim search for the oposite primer in each read, \
                    you should ensure it is never found at any other location than is intended. This is most \
                    problematic when your target can be very short and the oposite primer is found within the \
                    Illumina sequencing primer. It is recommended that you input >8 bases of each primer. the \
                    sequence should end with the last 3' base of the primer.")

        #Remove a primer from the primer set list
        elif sys.argv[1].upper() == 'REMOVE':
            if len(sys.argv) == 3:
                RemovePrimer(sys.argv[2].upper())
                exit()
            else:
                print("To remove a primer set from the list input the following argument variables: \
                    1. 'Remove'. \
                    2. Primer_Set_Name")

        elif len(sys.argv) != 8:
            print("To use input the following argument variables:\n\
                1. primer set name or 'other' if it is not in the common primer set list. \
                If a primer set name is entered, argument variables 2 & 3 will be ignored, but a value must be entered.\n\
                2. Last N bases of forward primer (recomended >= 8)\n\
                3. Last N bases of reverse primer (recomended >= 8)\n\
                4. N errors allowed in forward primer (recomended >=1)\n\
                5. N errors allowed in reverse primer (recomended >=1)\n\
                6. Length of marker. If you wish to use the length associated with a primer set, input the primer \
                set name. If length is variable, enter 0 and MetaTrim will search for the opposite primer in each \
                forward and reverse read. If the primer is not found it will take the remainder of the sequence after \
                the first primer is found. WARNING: if your marker is much shorter than the read length, you should \
                ensure that the oposite primer sequence is not found in the Illumina sequencing primer or in any \
                region after it's intended location.\n7. Are you using spacer inserts as published in Klymus et al. \
                2018, Plos One?: Y/N\n\n\
                To add a new primer set to the primer set list input the following argument vaiables: \
                1. 'New'. \
                2. Primer_Set Name, \
                3. Forward sequence, \
                4. Reverse sequence, \
                5. Length of marker OR primer set name to use a predefined length OR if length is variable, \
                enter 0 and MetaTrim will search for the opposite primer in each forward and reverse read. \
                If the primer is not found it will take the remainder of the sequence after the first primer is found. \
                WARNING: if you choose to have MetaTrim search for the oposite primer in each read, you should \
                ensure it is never found at any other location than is intended. This is most problematic when \
                your target can be very short and the oposite primer is found within the Illumina sequencing primer. \
                It is recommended that you input >8 bases of each primer. the sequence should end with the last 3' \
                base of the primer.\n\
                To remove a primer set from the list input the following argument variables: \
                1. 'Remove'. \
                2. Primer_Set_Name")
            exit()
        else:
            if sys.argv[6].isnumeric():
                LengthMarker = int(sys.argv[6])
            else:
                LengthMarker = sys.argv[6].upper()
                if LengthMarker not in PrimerSets:
                    print("Length ", LengthMarker," is not in the primer sets list!")
                    exit()
        cwd = os.getcwd()
    # def MetaTrimAll(Dir, PrimerSetAll, PFAll, PRAll, ErrFAll, ErrRAll, TargetLenAll, SpacersAll)
        #Create output directory and summary files
        
        basenm = os.path.basename(cwd)
        outsumname = basenm+'TrimSummary.txt'
        outsum = open(outsumname, "w")
        if sys.argv[7].upper() == 'Y':
            outsum.write('Sample\tReads\tF Seqs Trimmed\tR Seqs Trimmed\tSeqs F & R Trimmed\tShort Seqs\t\
                Seqs w/ Correct Spacer Combo\tSeqs w/o Correct Spacer Combo\n')
            outsum.close()
        else:
            outsum.write('Sample\tReads\tF Seqs Trimmed\tR Seqs Trimmed\tSeqs F & R Trimmed\tShort Seqs\n')
            outsum.close()
        
        #Create results directory
        ResDirName = basenm+'TrimmedFastqs'
        try:
            os.mkdir(ResDirName)
            print("Directory " , ResDirName ,  " Created ") 
        except FileExistsError:
            pass

        #Read in items in parent directory
        ItemsInParentDir = os.listdir(".")

        for x in ItemsInParentDir:
            #Open all subdirectories in parent directory
            if os.path.isdir(x):
                ItemsInSubDirs = sorted(os.listdir(x))
                #Find fastq files
                for y in ItemsInSubDirs:
                    if re.search('fastq\.gz$', y):
                        sample = x[0:re.search('_', y).start()]
                        #Find forward fastq file
                        if re.search('R1_001\.fastq\.gz', y):
                            InF = x+'/'+y
                        elif re.search('R2_001\.fastq\.gz', y):
                            InR = x+'/'+y
                            MetaTrim(InF, InR, sys.argv[1].upper(), sys.argv[3].upper(), sys.argv[4].upper(), \
                                sys.argv[4], sys.argv[5], LengthMarker, sys.argv[7].upper())
        outsum.close()

        end = datetime.now().time()
        print ('MetaTrim start:', start, '\n', 'MetaTrim end:', end)
    
