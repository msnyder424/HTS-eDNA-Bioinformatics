####MetaTrim####

#by M. R. Snyder 2018. 
#Written for Python 3
#MetaTrim trims paired end metabarcoding reads returned from Illumina sequencing for subsequent merging in other programs such as Dada2, Unoise, or OBITools. MetaTrim works on the parent directory level. It should be run from the parent directory in which all of your subdirectories with Illumina sequencing results are stored. MetaTrim handles zipped fastq.gz files. It will deposit your trimmed meta-barcoding results as unzipped fastq files into a new subdirectory named with your parent directory name and 'TrimmedFastqs'. If you have used primers with identical spacer inserts to those published in Klymus et al. 2017 in Plos One, you can use those inserts to remove instances of index hops. The most effective method of library prep for removing index hops is to combine every forward spacer with every reverse spacer.
#To use MetaTrim input the following argument variables: 
#1. primer set name or "other" if it is not in the common primer set list. If a primer set name is entered, argument variables 2 & 3 will be ignored, but a value must be entered.
#2. Last N bases of forward primer 
#3. Last N bases of reverse primer
#4. N errors allowed in forward primer
#5. N errors allowed in reverse primer 
#6. Length of marker OR primer set name to use a predefined length OR if length is variable, enter 0 and MetaTrim will search for the opposite primer in each forward and reverse read. If the primer is not found it will take the remainder of the sequence after the first primer is found. WARNING: if you choose to have MetaTrim search for the opposite primer in each read, you should ensure it is never found at any other location than is intended.
#7. 'Are you using spacer inserts as published in Klymus et al. 2018, Plos One?: Y/N'
#To add a new primer set to the primer set list input the following argument vaiables: 1. 'New' in primer set. 2. Primer_Set Name, 3. forward sequence, 4. reverse sequence, 5. Length of the marker or input 0 (zero) if you want MetaTrim to always search for the opposite primer in each read. It is recommended that you input >8 bases of each primer. the sequence should end with the last 3' base of the primer.
#To remove a primer set from the list input the following argument variables: 1. Remove. 2. Primer_Set_Name

import re
import numpy as np
from itertools import permutations
import sys
import gzip
import os
import operator
from datetime import datetime

start = datetime.now().time()

#BeginPrimerSetList: this list contains commonly used primer sets. Each key is the primer set name. The value has structure 'ForwardSequence-ReverseSequence-MarkerLength' See MetaTrimREADME for more info.

PrimerSets = {
    'GOBYPART': 'TWAAAATYGC-ACRTCWCGRC-167',
    'CARPPART': 'CYCTHCTAGG-CYCCRTTRGC-136',
    'CARPCOMPLETE': 'TGATGAAAYTTYGGMTCYCTHCTAGG-AARAAGAATGATGCYCCRTTRGC-136',
    'GOBYCOMPLETE': 'AACVCAYCCVCTVCTWAAAATYGC-AGYCANCCRAARTTWACRTCWCGRC-165',
    }

#EndPrimerSetList

IUPACAmb = {'R' : '[AG]', 'Y' : '[CT]', 'S' : '[GC]', 'W' : '[AT]', 'K' : '[GT]', 'M' : '[AC]', 'B' : '[CGT]', 'D' : '[AGT]', 'H' : '[ACT]', 'V' : '[ACG]', 'N' : '[ATCG]'}

inputs = sys.argv

try:
    PrimerSet = sys.argv[1].upper()
except IndexError:
    pass
try:
    if sys.argv[6].isnumeric():
        TargetLen = int(sys.argv[6])
        if TargetLen == 0:
            def RevComp(Seq):
                global SeqRC
                trans = str.maketrans('ATGC', 'TACG')
                SeqRC = Seq[::-1].translate(trans)
    else:
        TargetLen = sys.argv[6].upper()        
except IndexError:
    pass

if PrimerSet == 'NEW':
    NewOrNot = input("Would you like to rename this version of MetaTrim?(Y/N)").upper()
    if NewOrNot == 'Y':
        NewPyFileName = input("Name the new MetaTrim program file. Don't forget .py!")
    else:
        NewPyFileName = sys.argv[0]
    if sys.argv[2].upper() in PrimerSets:
        print (sys.argv[2].upper(), 'is already in the primer set list!')
        exit()
    else:
        print ('Adding primer set:', sys.argv[2].upper(), 'F seq:', sys.argv[3], 'R seq:', sys.argv[4], 'Target length:', sys.argv[5])
        with open('MetaTrim.py', 'r') as pyfileOld:
            lines = pyfileOld.readlines()
        for indx, val in enumerate(lines):
            if re.search('#BeginPrimerSetList', val):
                NewPLineIndx = indx + 3
                lines.insert(NewPLineIndx, "\t'"+sys.argv[2].upper()+"': '"+sys.argv[3]+"-"+sys.argv[4]+"-"+sys.argv[5]+"',\n")
                break
        with open(NewPyFileName, 'w') as pyfileNew:
            for i in lines:
                pyfileNew.write(i)
        exit()
elif PrimerSet == 'REMOVE':
    NewOrNot = input("Would you like to rename this version of MetaTrim?(Y/N)").upper()
    if NewOrNot == 'Y':
        NewPyFileName = input("Name the new MetaTrim program file. Don't forget .py!")
    else:
        NewPyFileName = sys.argv[0]
    if not sys.argv[2].upper() in PrimerSets:
        print (sys.argv[2].upper(), 'is not in the primer set list!')
        exit()
    else:
        PrimerToRemove = "\t'"+sys.argv[2].upper()+"':"
        print ('Removing primer set', sys.argv[2].upper(), PrimerToRemove)
        with open('MetaTrim.py', 'r') as pyfileOld:
            lines = pyfileOld.readlines()
            for i in lines:
                if re.match(PrimerToRemove, i):
                    lines.remove(i)
                if re.match('#EndPrimerSetList', i):
                    break
        with open(NewPyFileName, 'w') as pyfileNew:
            for i in lines:
                pyfileNew.write(i)
        exit()
elif len(inputs) != 8:
    print("to use input the following argument variables:\n1. primer set name or 'other' if it is not in the common primer set list. If a primer set name is entered, argument variables 2 & 3 will be ignored, but a value must be entered.\n2. Last N bases of forward primer (recomended >= 8)\n3. Last N bases of reverse primer (recomended >= 8)\n4. N errors allowed in forward primer (recomended >=1)\n5. N errors allowed in reverse primer (recomended >=1)\n6. Length of marker. If you wish to use the length associated with a primer set, input the primer set name. If length is variable, enter 0 and MetaTrim will search for the opposite primer in each forward and reverse read. If the primer is not found it will take the remainder of the sequence after the first primer is found. WARNING: if your marker is much shorter than the read length, you should ensure that the oposite primer sequence is not found in the Illumina sequencing primer or in any region after it's intended location.\n7. Are you using spacer inserts as published in Klymus et al. 2018, Plos One?: Y/N\n\nTo add a new primer set to the primer set list input the following argument vaiables: 1. 'New' in primer set. 2. Primer_Set Name, 3. Forward sequence, 4. Reverse sequence, 5. Length of marker OR primer set name to use a predefined length OR if length is variable, enter 0 and MetaTrim will search for the opposite primer in each forward and reverse read. If the primer is not found it will take the remainder of the sequence after the first primer is found. WARNING: if you choose to have MetaTrim search for the oposite primer in each read, you should ensure it is never found at any other location than is intended. This is most problematic when your target can be very short and the oposite primer is found within the Illumina sequencing primer. It is recommended that you input >8 bases of each primer. the sequence should end with the last 3' base of the primer.\nTo remove a primer set from the list input the following argument variables: 1. Remove. 2. Primer_Set_Name")
    exit()
elif PrimerSet == 'OTHER':
    print ('Primer set is ', PrimerSet,'. F: ', sys.argv[2].upper(), 'R: ', sys.argv[3].upper())
    PrimerDict = {'pF': sys.argv[2].upper(), 'pR': sys.argv[3].upper()}
    ErrDict = {'pF': int(sys.argv[4]), 'pR': int(sys.argv[5])}
    TargetLen = int(sys.argv[6])
elif PrimerSet in PrimerSets:
    print ('Primer set is ', PrimerSet)
    Vars = PrimerSets[PrimerSet].split('-')
    PrimerDict = {'pF': Vars[0], 'pR': Vars[1]}
    ErrDict = {'pF': int(sys.argv[4]), 'pR': int(sys.argv[5])}
    if TargetLen in PrimerSets:
        TargetLen = int(Vars[2])
    elif TargetLen > 0:
        TargetLen = int(sys.argv[6])
elif PrimerSet not in PrimerSets:
    print (PrimerSet,'is not in the Primer Sets List!')
    exit()

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
        if TargetLen == 0:
            RevComp(readbuffer[1])
            if re.search(DegPrimerDict[Primer2], SeqRC):
                TargetEnd = len(readbuffer[1]) - re.search(DegPrimerDict[Primer2], SeqRC).end()
            else:
                TargetEnd = len(readbuffer[1])
        else:
            TargetEnd = TargetStart + TargetLen
        Seq = readbuffer[1][TargetStart:TargetEnd]
        Qual = readbuffer[3][TargetStart:TargetEnd]

if sys.argv[7] == 'Y':
    FSpacers = {'e': '(TCCT|[ATCGN]CCT|T[ATCGN]CT|TC[ATCGN]T|TCC[ATCGN])', 'f': '(ATGC|[ATCGN]TGC|A[ATCGN]GC|AT[ATCGN]C|ATG[ATCGN])', 'g': '(CGAG|[ATCGN]GAG|C[ATCGN]AG|CG[ATCGN]G|CGA[ATCGN])', 'h': '(GATA|[ATCGN]ATA|G[ATCGN]TA|GA[ATCGN]A|GAT[ATCGN])'}
    RSpacers = {'e': '(CGTA|[ATCGN]GTA|C[ATCGN]TA|CG[ATCGN]A|CGT[ATCGN])', 'f': '(TCAC|[ATCGN]CAC|T[ATCGN]AC|TC[ATCGN]C|TCA[ATCGN])', 'g': '(GAGT|[ATCGN]AGT|G[ATCGN]GT|GA[ATCGN]T|GAG[ATCGN])', 'h': '(ATCG|[ATCGN]TCG|A[ATCGN]CG|AT[ATCGN]G|ATC[ATCGN])'}
    def SpacerCount(Spacer, ReadDirection, Dir):
        Count = {'e': 0, 'f': 0, 'g': 0, 'h': 0}
        print ('Removing index hops based on incorrect spacer insert in sample', Dir, ReadDirection)
        readbuffer = []
        global reads
        global MaxS
        reads = 0
        for line in infile:
            LineCurr = line.strip().decode('ascii')
            readbuffer.append(LineCurr)
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
        print (Dir, ReadDirection, 'has spacer', MaxS)

DegPrimerDict = {}
print ('Making degenerate primer regex.')
for keys in PrimerDict:
    DegPrimers = []
    PrimCurr = list(PrimerDict[keys])
    PrimDegDict = {}
    for err in range(0, ErrDict[keys]+1):
        PrimDeg = []
        for bp in range(len(PrimerDict[keys])):
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

cwd = os.getcwd()
basenm = os.path.basename(cwd)
outsumname = basenm+'TrimSummary.txt'
outsum = open(outsumname, "w")
if sys.argv[7] == 'Y':
    outsum.write('Sample\tReads\tF Seqs Trimmed\tR Seqs Trimmed\tSeqs F & R Trimmed\tShort Seqs\tSeqs w/ Correct Spacer Combo\tSeqs w/o Correct Spaer Combo\n')
else:
    outsum.write('Sample\tReads\tF Seqs Trimmed\tR Seqs Trimmed\tSeqs F & R Trimmed\tShort Seqs\n')
SubDirs = os.listdir(".")
ResDirName = basenm+'TrimmedFastqs'

try:
    os.mkdir(ResDirName)
    print("Directory " , ResDirName ,  " Created ") 
except FileExistsError:
    pass

for x in SubDirs:
    if os.path.isdir(x):
        sample = x[0:6]
        FNames = {}
        FSeqs = {}
        FQuals = {}
        CorSpacer = {}
        RSeqs = 0
        FinalSeqs = 0
        Shorts = 0
        files = sorted(os.listdir(x))
        for y in files:
            if re.search('fastq\.gz$', y):
                InF = x+'/'+y
                if re.search('R1_001\.fastq\.gz', y):
                    reads = 0
                    Foutfilenameandpath = ResDirName+'/'+sample+'_R1_001.fastq'
                    Foutfile = open(Foutfilenameandpath, "w")
                    readbuffer = []
                    print ('Processing:', InF)
                    with gzip.open(InF) as infile:
                        if sys.argv[7] == 'Y':
                            SpacerCount(FSpacers, 'forward', sample)
                            FMax = MaxS
                            print ('Trimming primers sample', x, 'forward.')
                            infile.seek(0, 0)
                            for line in infile:
                                LineCurr = line.strip().decode('ascii')
                                readbuffer.append(LineCurr)
                                if len(readbuffer) == 4:
                                    reads += 1
                                    if reads % 10000 == 0:
                                        print ('Read:', reads, end='\r')
                                    if re.match(FSpacers[FMax], readbuffer[1]):
                                        TrimPrimers('pF','pR')
                                        if len(Seq) > 100:
                                            FNames[name[0]] = readbuffer[0]
                                            FSeqs[name[0]] = Seq
                                            FQuals[name[0]] = Qual
                                            CorSpacer[name[0]] = 1
                                    readbuffer = []
                            print(reads, 'total raw reads.')
                        else:
                            print ('Trimming primers sample', x, 'forward.')
                            for line in infile:
                                LineCurr = line.strip().decode('ascii')
                                readbuffer.append(LineCurr)
                                if len(readbuffer) == 4:
                                    reads += 1
                                    if reads % 10000 == 0:
                                        print ('Read:', reads, end='\r')
                                    TrimPrimers('pF','pR')
                                    if len(Seq) > 100:
                                        FNames[name[0]] = readbuffer[0]
                                        FSeqs[name[0]] = Seq
                                        FQuals[name[0]] = Qual
                                    readbuffer = []
                            print(reads, 'total raw reads.')
                elif re.search('R2_001', y):
                    Routfilenameandpath = ResDirName+'/'+sample+'_R2_001.fastq'
                    Routfile = open(Routfilenameandpath, "w")
                    readbuffer = []
                    reads = 0
                    print ('opening:', InF)
                    with gzip.open(InF) as infile:
                        if sys.argv[7] == 'Y':
                            SpacerCount(RSpacers, 'reverse', sample)
                            RMax = MaxS
                            print ('Trimming primers sample', x, 'reverse.')
                            infile.seek(0, 0)
                            for line in infile:
                                LineCurr = line.strip().decode('ascii')
                                readbuffer.append(LineCurr)
                                if len(readbuffer) == 4:
                                    reads += 1
                                    if reads % 10000 == 0:
                                        print ('Read:', reads, end='\r')
                                    if re.match(RSpacers[RMax], readbuffer[1]):
                                        TrimPrimers('pR','pF')
                                        RSeqs += 1
                                        if name[0] in FSeqs:
                                            if len(Seq) > 100:
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
                            print ('Trimming primers sample', x, 'reverse.')
                            for line in infile:
                                LineCurr = line.strip().decode('ascii')
                                readbuffer.append(LineCurr)
                                if len(readbuffer) == 4:
                                    reads += 1
                                    if reads % 10000 == 0:
                                        print ('Read:', reads, end='\r')
                                    TrimPrimers('pR','pF')
                                    RSeqs += 1
                                    if name[0] in FSeqs:
                                        if len(Seq) > 100:
                                            Foutfile.write("%s\n%s\n+\n%s" % (FNames[name[0]],FSeqs[name[0]],FQuals[name[0]]))
                                            Routfile.write("%s\n%s\n+\n%s" % (readbuffer[0],Seq,Qual))
                                            FinalSeqs += 1
                                        else:
                                            Shorts += 1
                                    readbuffer = []
                            print(reads, 'total raw reads.')
                    Foutfile.close()
                    Routfile.close()
                    if sys.argv[7] == 'Y':
                        CorrSpacerCount = sum(CorSpacer.values())
                        IncorrSpacerCount = reads - CorrSpacerCount
                        outsum.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sample, reads, len(FSeqs), RSeqs, FinalSeqs, Shorts, CorrSpacerCount, IncorrSpacerCount))
                    else:
                        outsum.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (sample, reads, len(FSeqs), RSeqs, FinalSeqs, Shorts))

end = datetime.now().time()
print ('MetaTrim start:', start, '\n', 'MetaTrim end:', end)
