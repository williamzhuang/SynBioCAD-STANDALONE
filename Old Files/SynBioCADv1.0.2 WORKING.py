#!/usr/bin/env python
#title           : SynBioCAD1.0 beta
#description     : defines several methods to generate oligos for CRISPR
#author          : Leo d'Espaux <leodespaux@gmail.com>
#date            : 9 June 2015
#==============================================================================


# import libraries we're using
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp
import copy

from intermine.webservice import Service

from pandas import *
from pandas import DataFrame, read_csv
import pandas as pd #this is how I usually import pandas

import matplotlib.pyplot as plt

import sys

# define global variables
HomologyLength = 1000
PrimerMaxTm = 55
PrimerMaxLen = 60
OverhangMaxFrac = 1





def askUser():
    print('Python version ' + sys.version)
    print("IMPORTANT: IF THE ABOVE SAYS PYTHON 3.X DO NOT USE QUOTES IN RESPONSES.")
    #print('Pandas version ' + pd.__version__)
    print(" ")
    print("Hi baby, welcome to SynBioCAD. I'm here for you.")
    print("Right now, we assume you're working in S. cerevisiae.")
    print(" ")
    print("What do you want to do today?")
    print("\"1\" to insert DNA into a characterized locus") 
    print("\"2\" to edit an existing gene, e.g., delete or replace")
    print("\"3\" to build a custom cassette")
    print(" ")

    Action= input("Your answer: ")
    if Action == "1":
        print(" ")
        editEmpty()
        
    elif Action == "2":
        print(" ")
        editExisting()

    elif Action == "3":
        print(" ")
        buildCustom()
              
def editEmpty():
    
    labels=['208a', '1014a', '1114a', '607c', '308a', '1021b', '720a']
    ChrLetters=["Scer02", "Scer10", "Scer11", "Scer06", "Scer03", "Scer10", "Scer07"]
    ExpValues=[1.0, 1.2, 1.1, 1.1, 1.4, 1.0, 1.0]
    cutSeqs=["GTCCGCTAAACAAAAGATCT", "TTATGTGCGTATTGCTTTCA", "CTTGTGAAACAAATAATTGG", "CTATTTTTGCTTTCTGCACA", "TAGGATACAAGTGTCTGTGC", "CCTCTGTGTGGTGGTAATTG", "CAACAATTGTTACAATAGTA"]

    
    #We save information as an array
    cutArray={'name' : Series(labels, index=labels),
                'exp. lev.' : Series(ExpValues, index=labels),
                'chrom. loc.' : Series(ChrLetters, index=labels),
                'sequence' : Series(cutSeqs, index=labels)
             }

    
    cutFrame= DataFrame(cutArray)
    
    print(cutFrame)
    
    print("")
    print("Which cut site do you want, e.g., 208a?")
    cutname=input("Your answer: ")

    location=cutFrame.loc[cutname,'chrom. loc.']+".fasta"
    cutSequence=cutFrame.loc[cutname,'sequence']
    
    
    # we use that location to get the chromosome sequence
    ChromosomeSeq=SeqIO.read(location, "fasta").seq
    
    if ChromosomeSeq.find(cutSequence)==-1:
        ChromosomeSeq=ChromosomeSeq.reverse_complement()

    StartIndex=ChromosomeSeq.find(cutSequence)
    EndIndex=StartIndex+34
    
    UpSeq=ChromosomeSeq[StartIndex-HomologyLength:StartIndex]
    DownSeq=ChromosomeSeq[EndIndex:EndIndex+HomologyLength]
        
    UpHomRec = SeqRecord(UpSeq, id=cutname)
    DownHomRec = SeqRecord(DownSeq, id=cutname)
    
    print("")
    print("Found your homology regions")
    print("")

              
    print("What do you want to do to " + cutname + "?")
    print("\"1\" to insert a pre-built custom cassette")
    print("\"2\" to construct a cassette using standard promoter and terminator names, and a custom ORF")
    print(" ")

    typeEdit=input("Your answer: ")
    if typeEdit=="1":
    
        orfName=input("What's the name of your custom gene or cassette?")
        orfSeq=Seq(input("What's the sequence?"))
        
        orfRecord=SeqRecord(orfSeq,id=orfName)
        
        fragments=[UpHomRec, orfRecord, DownHomRec]

    
    elif typeEdit == "2":
        
        PromoterRec, orfRecord, TerminatorRec = buildCassette()
        
        fragments=[UpHomRec, PromoterRec, orfRecord, TerminatorRec, DownHomRec]

    
    stitch(fragments)
    
def editExisting():
    print(" ")
    print("Which locus do you want to edit? Tell me a common name, e.g., \"OAF1\": ")
    print("I'm smart and pretty and I can fetch it for you.")
    
    GeneName= input("Your answer: ")
    
    print("")
    print("Give me a second to find it for you...")

        
    OrigGeneRecord=fetchGene(GeneName)
    #note that this returns a seqrecord
    
    print("")
    print("I think you want "+OrigGeneRecord.features+":")
    print(OrigGeneRecord.description)


    # We make seqrecords since that's what we carry through later in the program
    UpHomRec = fetchNeighbor(OrigGeneRecord, "upstream", HomologyLength )
    DownHomRec = fetchNeighbor(OrigGeneRecord, "downstream", HomologyLength )
    
    #print("OK I found them.")

    #print(" ")
    #print("UpstreamHomology reghhion is")
    #print(UpHomRec.seq)
    
    #print(" ")
    #print("DownstreamHomology reghhion is")
    #print(DownHomRec.seq)
    
    print(" ")
    print("What do you want to do to this gene?")
    print("\"1\" delete the protein coding DNA sequence (CDS)")
    print("\"2\" replace the CDS with another CDS, or other pre-made fragment")
    print("\"3\" replace the CDS with a standard cassette I will build for you")
    print("\"4\" replace the CDS with a custom cassette")
    print("\"5\" replace a specified region near your target gene")
    
    print(" ")
    Action=input("Your answer: ")
                        
    #note that in all the below, we want to have fragments be records
    if Action=="1":
        fragments=[UpHomRec,DownHomRec]
    
    if Action=="2":
        print(" ")
        NewGeneName=input("What's the name of the gene you're inserting?")
        NewGeneSeq=Seq(input("What's the sequence of your new gene?"))
        InsertRec = SeqRecord(NewGeneSeq, id=NewGeneName)
        fragments=[UpHomRec, InsertRec, DownHomRec]

    if Action=="3":
        print(" ")
        print("I got you baby. I don't want you to strain your eyes looking up stuff.")
        PromoterRec, orfRecord, TerminatorRec = standardCassette()
        fragments=[UpHomRec, PromoterRec, orfRecord, TerminatorRec, DownHomRec] #we need to finish buildcassette to add InsertRec here

    if Action=="4":
        print("How many pieces (other than homology fragments) are you stitching together.")
        Npieces=input("Your answer: ")
        output = standardCassette()

        fragments=[UpHomRec, DownHomRec] #we need to finish buildcassette to add InsertRec here

        
    stitch(fragments)

def buildCustom():
    print("How do you want your cassettes?")
    print("\"1\" to build a cassette with N pieces and no variants.")
    print("\"2\" to build a cassette with N pieces and one variant site.")

    answer = input("Your answer:")

    if answer == "1":
        N = int(input("How many pieces (N) in your custom cassette:"))
        fragments = variableCassette(N)[0]
        stitch(fragments)


    elif answer == "2":
        N = int(input("How many pieces (N) are in your custom cassette: "))
        toVary = int(input("Which piece would you like to vary: "))
        variants = int(input("How many variants would you like: "))
        fragsList = variableCassette(N, toVary, variants)
        # Make the entire thing first. Then just primers for the variant.
        print("")
        print("First we'll get primers for the entire cassette")
        stitch(fragsList[0])
        print("-------------------------------------------------")
        print("Now we'll get primers just for the variants")

        for n in range(N-1):
            print("For variant " + str(n) + ":")
            frags = fragsList[n + 1][toVary - 2: toVary + 1]
            stitch(frags)
    #------------------------------ FETCH FUNCTIONS -------------------------------------

def fetchGene(GeneName):
    
    service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")
    template = service.get_template('Gene_GenomicDNA')

    rows = template.rows(
        E = {"op": "LOOKUP", "value": GeneName, "extra_value": "S. cerevisiae"}
    )
    
    # this service seems to return multiple similar genes but we want the first one only, so count
    # and it returns information about the gene you want
    count=0
    for row in rows:
        
        count=count+1
        if count==1:
            descr= row["description"]
            GeneSeq=Seq(row["sequence.residues"])
            GeneSysName=row["secondaryIdentifier"]
            #print(" ")
            #print("I think you want...... "+row["secondaryIdentifier"])
            #print(row["description"])
            #print(" ")
            #print(row["sequence.residues"])
            #print(" ")
            #print("Good choice! I have a feeling you're going to get lucky with this one.")
            #print(" ")
            #print("Give me a second to put some of my ducks in a circle...")
       

            
    #let's create a record for the oldGene
    GeneRecord = SeqRecord(GeneSeq, id=GeneSysName)
    
    #now let's add some more information to make it useful
    GeneRecord.name=GeneName
    GeneRecord.features=GeneSysName
    GeneRecord.description=descr

    return GeneRecord 
       
def fetchNeighbor(NeighborRecord, direction, distance):


    # let's load the appropriate chromosome file. The record of the gene we looked up
    # contains in the "features" the systematic name, wherein the second letter
    # corresponds to chromosome number, e.g., 1=A etc
    if NeighborRecord.features[1]=="A":
        ChromosomeRec=SeqIO.read("Scer01.fasta", "fasta")
    if NeighborRecord.features[1]=="B":
        ChromosomeRec=SeqIO.read("Scer02.fasta", "fasta")
    if NeighborRecord.features[1]=="C":
        ChromosomeRec=SeqIO.read("Scer03.fasta", "fasta")
    if NeighborRecord.features[1]=="D":
        ChromosomeRec=SeqIO.read("Scer04.fasta", "fasta")
    if NeighborRecord.features[1]=="E":
        ChromosomeRec=SeqIO.read("Scer05.fasta", "fasta")
    if NeighborRecord.features[1]=="F":
        ChromosomeRec=SeqIO.read("Scer06.fasta", "fasta")
    if NeighborRecord.features[1]=="G":
        ChromosomeRec=SeqIO.read("Scer07.fasta", "fasta")
    if NeighborRecord.features[1]=="H":
        ChromosomeRec=SeqIO.read("Scer08.fasta", "fasta")
    if NeighborRecord.features[1]=="I":
        ChromosomeRec=SeqIO.read("Scer09.fasta", "fasta")
    if NeighborRecord.features[1]=="J":
        ChromosomeRec=SeqIO.read("Scer10.fasta", "fasta")
    if NeighborRecord.features[1]=="K":
        ChromosomeRec=SeqIO.read("Scer11.fasta", "fasta")
    if NeighborRecord.features[1]=="L":
        ChromosomeRec=SeqIO.read("Scer12.fasta", "fasta")
    if NeighborRecord.features[1]=="M":
        ChromosomeRec=SeqIO.read("Scer13.fasta", "fasta")
    if NeighborRecord.features[1]=="N":
        ChromosomeRec=SeqIO.read("Scer14.fasta", "fasta")
    if NeighborRecord.features[1]=="O":
        ChromosomeRec=SeqIO.read("Scer15.fasta", "fasta")
    if NeighborRecord.features[1]=="P":
        ChromosomeRec=SeqIO.read("Scer16.fasta", "fasta") 

    
    
    # let's explicitely name the sequences from the seq record
    NeighborSeq=NeighborRecord.seq
    ChromosomeSeq=ChromosomeRec.seq
    
    # flip the sequence to orient with respect to the old gene
    if ChromosomeSeq.find(NeighborSeq)==-1:
        ChromosomeSeq=ChromosomeSeq.reverse_complement()

    StartIndex=ChromosomeSeq.find(NeighborSeq)
    EndIndex=StartIndex+len(NeighborSeq)
    
    if direction=="upstream":
        DesiredSeq=ChromosomeSeq[StartIndex-distance:StartIndex]
    if direction=="downstream":
        DesiredSeq=ChromosomeSeq[EndIndex:EndIndex+distance]

    
    
    
    NeighborRec = SeqRecord(DesiredSeq, id=NeighborRecord.name)
    
    return NeighborRec

    #print(NeighborRec)

    
    
    
    
    
    
    
    
    
    
    
    #------------------------------------ CONSTRUCTING STUFF --------------------------------------

def getPrimer(currRecord):
    

    mp = 0
    length = 0
    primer = Seq("")

    seq=currRecord.seq
    
    while mp <= PrimerMaxTm and length <= PrimerMaxLen:
        primer = primer + seq[length]
        mp = MeltingTemp.Tm_staluc(primer)
        length += 1

    return primer           
        
def overhangPrimer(currRecord,prevSeq):
    #let's get the template-binding primer first
    primer=getPrimer(currRecord)
    
    
    #OK let's work on the overhang
    maxOhLen=PrimerMaxLen-len(primer)    
    maxFrac=1
    
    #let's decide on a max overhang length
    if round(len(primer)*(OverhangMaxFrac+1)) < 60:
             maxOhLen=round(len(primer)*OverhangMaxFrac)
    
    #the index must be an integer!!!
    maxOhLen=int(maxOhLen)
    ohprimer=prevSeq.seq[-maxOhLen:]+primer #we add the .seq so that it returns a string
    
    return ohprimer      
        
def standardCassette():
    
    #first, the promoter
    print("I'm going to build a standard cassette in which promoter is 600nt, terminator 250nt.") 
    print("First, which PROMOTER do you want to use, e.g., TDH3")
    
    PromoterName=input("Your answer: ")
    PromoterGeneRec=fetchGene(PromoterName)
    PromoterRec=fetchNeighbor(PromoterGeneRec,"upstream",600)
    PromoterRec.id=PromoterRec.id+"ps"
    
    
    #second, the terminator
    print("Which TERMINATOR do you want to use, e.g., ADH1")
    TerminatorName = input('Your answer: ')
    TerminatorGeneRec=fetchGene(TerminatorName)
    TerminatorRec=fetchNeighbor(TerminatorGeneRec,"downstream",250)
    TerminatorRec.id=TerminatorRec.id+"ts"
    
    
    #and last, the gene
    print("What is the name of your gene, e.g., KlGapDH")
    orfName = input("Your answer: ")
    
    print("What's the sequence")
    orfSeq=input("Your answer: ")
    
    orfRecord=SeqRecord(Seq(orfSeq), id=orfName)
    
    insertRec=[PromoterRec,orfRecord,TerminatorRec]
    return PromoterRec, orfRecord, TerminatorRec
          
def buildCassette():
    
    #first, the promoter
    print("I'm going to build a standard cassette in which promoter is 600nt, terminator 250nt.") 
    print("First, which PROMOTER do you want to use, e.g., TDH3")
    
    PromoterName=input("Your answer: ")
    PromoterGeneRec=fetchGene(PromoterName)
    PromoterRec=fetchNeighbor(PromoterGeneRec,"upstream",600)
    PromoterRec.id=PromoterRec.id+"ps"
    
    
    #second, the terminator
    print("Which TERMINATOR do you want to use, e.g., ADH1")
    TerminatorName = input('Your answer: ')
    TerminatorGeneRec=fetchGene(TerminatorName)
    TerminatorRec=fetchNeighbor(TerminatorGeneRec,"downstream",250)
    TerminatorRec.id=TerminatorRec.id+"ts"
    
    #and last, the gene
    print("What is the name of your gene, e.g., KlGapDH")
    orfName = input("Your answer: ")
    
    print("What's the sequence")
    orfSeq=input("Your answer: ")
    
    orfRecord=SeqRecord(Seq(orfSeq), id=orfName)
    
    insertRec=[PromoterRec,orfRecord,TerminatorRec]
    return PromoterRec, orfRecord, TerminatorRec
    
def variableCassette(N, toVary = 0, variants = 0):
    print("")
    print("Let's start building.")
    print("")
    
    # Store both name and sequence in a SeqRecord
    # Append them to a list
    # Return list as fragments to be stitched

    records = []
    for n in range(N):
        name = input("What is the name of sequence " + str(n+1) +":")
        sequence = input("What is the sequence of this fragment:")
        print("")
        Rec = SeqRecord(Seq(sequence), id = name)
        Rec.name = name
        records.append(Rec)

    variantRecords = []
    variantRecords.append(records)
    # This only happens if there are variants.
    if variants > 0:
        print("Time to make those variants you wanted.")
        for n in range(variants-1):
            name = input("What is the name of variant " + str(n+1) + ":")
            sequence = input("What is the sequence of this variant:")
            Rec = SeqRecord(Seq(sequence), id = name)
            Rec.name = name
            # Make a copy of the original, switch the fragments and add it to the list. 
            # Deep-copy ensures there are no pointer issues
            tempVariant = copy.deepcopy(records)
            tempVariant[toVary - 1] = Rec
            variantRecords.append(copy.deepcopy(tempVariant))
            print("")

    # Returns a list of lists of the SeqRecords of the fragments
    return variantRecords

def stitch(fragments):
    #this function takes seq records and prints primers
    
    #let's make an empty sequence file
    Nfrags=len(fragments)
    donor=Seq("")
    index=[]
    print("")
    for i in range (0, Nfrags):
        donor=donor+fragments[i]
    
    for i in range (0, Nfrags):
        if i==0:
            print("Lup"+ fragments[i].id + " " + getPrimer(donor))
            print("Rup"+ fragments[i].id + "(" + fragments[i+1].id + ") " + overhangPrimer(fragments[i].reverse_complement(),fragments[i+1].reverse_complement()))
        elif i==Nfrags-1:
            print("Ldown"+ fragments[i].id + "(" + fragments[i-1].id + ") " + overhangPrimer(fragments[i],fragments[i-1]))
            print("Rdown"+ fragments[i].id + " " + getPrimer(donor.reverse_complement()))
        else:
            print("L"+ fragments[i].id + "(" + fragments[i-1].id + ") " + overhangPrimer(fragments[i],fragments[i-1]))
            print("R"+ fragments[i].id + "(" + fragments[i+1].id + ") " + overhangPrimer(fragments[i].reverse_complement(),fragments[i+1].reverse_complement()))

    print("")
    print("Your donor DNA cassette, has the following bp length and sequence:")


    print("")
    print(len(donor.seq))
    print("")

    print(donor.seq)

    print("")
    print("You might want to copy this entire prompt and save it for your records.")
        

    
def pickCut(Sequence):
    from prettytable import from_csv

    print("Which of these characterized loci is going to get lucky: ")

    fp = open("cutsites.csv", "r")
    pt = from_csv(fp)
    fp.close()

    print(pt)

    cutsiteNum = int(input("Choose your cutsite (Leftmost number):"))

    print("")
    print("You've chosen: ")
    print(pt.get_string(start = cutsiteNum-1, end = cutsiteNum))

    counter = 0
    for row in pt:
        counter += 1
        if counter == cutsiteNum:
            row.border = False
            row.header = False
            cutsite = row.get_string(fields=["sequence"]).strip()
            cutName = row.get_string(fields=["name"]).strip()
            chromosome = int(row.get_string(fields=["chromosome"]).strip())
    if chromosome >= 10:
        filename = "Scer" + str(chromosome) + ".fasta" 
    else:
        filename = "Scer0" + str(chromosome) + ".fasta"

    ChromosomeRec=SeqIO.read(filename, "fasta")
    ChromosomeSeq = ChromosomeRec.seq

    if ChromosomeSeq.find(cutsite) == -1:
        ChromosomeSeq = ChromosomeSeq.reverse_complement()

    startInd = ChromosomeSeq.find(cutsite)

    UpHomology = SeqRecord(ChromosomeSeq[startInd - 1000:startInd])
    DownHomology = SeqRecord(ChromosomeSeq[startInd + 30:startInd + 1030])
    return UpHomology, DownHomology



#fetchGene("TEF1")   
#print (fetchNeighbor("TEF1", "downstream", 250).seq)

#print (buildCassette())
askUser()

