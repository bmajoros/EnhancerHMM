#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import ProgramName
#import commands
from Pipe import Pipe

def fisher(a,b,c,d):
    cmd="/home/bmajoros/src/scripts/fisher-exact-test.R "+str(a)+" "+\
        str(b)+" "+str(c)+" "+str(d)
    #(ret,line)=commands.getstatusoutput(cmd)
    line=Pipe.run(cmd)
    return float(line.rstrip())

def getTable(motif1,motif2,peaks,singletons):
    yesyes=0; yesno=0; noyes=0; nono=0
    for pair in peaks:
        (peakID,peak)=pair
        if(singletons.get(peakID,None) is None): continue
        if(pairPresent(motif1,motif2,peak)):
            yesyes+=1
            continue
        present=singletons[peakID]
        motif1present=motif1 in present
        motif2present=motif2 in present
        if(motif1==motif2):
            if(motif1present):
                yesno+=1; noyes+=1
            else: nono+=1
        else:
            if(motif1present): yesno+=1
            if(motif2present): noyes+=1
            if(not (motif1present or motif2present)): nono+=1
    return (yesyes,yesno,noyes,nono)

def getTable_old(motif1,motif2,peaks):
    yesyes=0; yesno=0; noyes=0; nono=0
    for peak in peaks:
        if(pairPresent(motif1,motif2,peak)):
            yesyes+=1
            continue
        present=getPresent(peak)
        motif1present=motif1 in present
        motif2present=motif2 in present
        if(motif1==motif2):
            if(motif1present):
                yesno+=1; noyes+=1
            else: nono+=1
        else:
            if(motif1present): yesno+=1
            if(motif2present): noyes+=1
            if(not (motif1present or motif2present)): nono+=1
    return (yesyes,yesno,noyes,nono)

def pairPresent(motif1,motif2,peak):
    for pair in peak:
        if(pair[0]==motif1 and pair[1]==motif2): return True
    return False

def getPresent(peak):
    present=set()
    for pair in peak:
        present.add(pair[0])
        present.add(pair[1])
    return present

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+
         " <motif-pairs.txt> <peaks-motifs.txt>\n")
(pairsFile,singletonsFile)=sys.argv[1:]

# Read singleton motif occurrences
singletons={} # peak -> array of motifs
IN=open(singletonsFile,"rt")
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=2): continue
    (peak,motif)=fields
    if(singletons.get(peak,None) is None): singletons[peak]=set()
    singletons[peak].add(motif)
IN.close()

# Read the file of motif pair occurrences
peaks=[]
pairs=[]
motifs=set()
peakID=None
IN=open(pairsFile,"rt")
IN.readline()
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=2): continue
    (first,second)=fields
    if(first=="peak"):
        peaks.append([peakID,pairs])
        pairs=[]
        peakID=second
    else: 
        pairs.append([first,second])
        motifs.add(first); motifs.add(second)
peaks.append([peakID,pairs])
IN.close()

# Compute expected counts and get P-values
for motif1 in motifs:
    for motif2 in motifs:
        if(motif1>motif2): continue
        (a,b,c,d)=getTable(motif1,motif2,peaks,singletons)
        P=fisher(a,b,c,d)
        #if(P<0.05):
        if(True):
            N=float(a+b+c+d)
            Ptop=float(a+b)/N
            Pbottom=1.0-Ptop
            expA=int(Ptop*float(a+c))
            expB=Ptop*float(b+d)
            expC=Pbottom*float(a+c)
            expD=Pbottom*float(b+d)
            #if(a<=expA): continue
            print(P,motif1,motif2,a,expA,sep="\t")
            print("\t"+str(a)+"\t"+str(b)+"\t"+str(c)+"\t"+str(d))
            print("\t"+str(expA)+"\t"+str(expB)+"\t"+str(expC)+"\t"+str(expD))
