#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import ProgramName
from Rex import Rex
rex=Rex()

SINGLE_PEAK_ONLY=False
MULTI_PEAK_ONLY=True

def loadLoopFile(filename):
    hash={} # maps genes to enhancers
    with open(filename) as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=2): continue
            (gene,peak)=fields
            if(hash.get(gene,None) is None): hash[gene]=[]
            hash[gene].append(peak)
    return hash

def loadP300file(filename,whichTime):
    # iter0_peak14413.t3.fastb        2234.2906200000034      1.0     1:0-1|2:1-2|3:2-1243|4:1243-1256|3:1256-1998|4:1998-1999|5:1999-2000    1.07942122577,1.12211048111,1.00477844983,1.16141454506|0.952633976593,1.02660980228,0.892990209161,0.961271510352      10      10      11      10      10      00
    hash={} # maps enhancerID to P300 value
    with open(filename) as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=11): continue
            (fastb,LLR,P,parse,features,GR,AP1,CEBP,FOX,KLF,CTCF)=fields
            if(not rex.find("(\S+)\.(t\d+)\.fastb",fastb)):
                raise Exception(fastb)
            enhancerID=rex[1]
            time=rex[2]
            if(time!=whichTime): continue
            fields=features.split("|")
            numPeaks=len(fields)
            if(SINGLE_PEAK_ONLY and numPeaks>1): continue
            if(MULTI_PEAK_ONLY and numPeaks<2): continue
            totalP300=0
            for field in fields:
                subfields=field.split(",")
                (DNase_t0,DNase_t3,P300_t0,P300_t3)=subfields
                delta=float(P300_t3)-float(P300_t0)
                totalP300+=delta
            hash[enhancerID]=totalP300
    return hash

def loadExpressionFile(filename):
    # gene        logFC   logCPM  LR      PValue  FDR
    hash={}
    with open(filename) as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=6): continue
            (gene,logFC,logCPM,LR,Pvalue,FDR)=fields
            hash[gene]=logFC
    return hash

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <genes-enhancers.txt> <genomewide-features.txt> <edgeR.txt> <enhancer-timepoint>\n")
(loopFile,p300file,expressionFile,enhancerTime)=sys.argv[1:]

geneToEnhancer=loadLoopFile(loopFile)
enhancerToP300=loadP300file(p300file,enhancerTime)
geneToFC=loadExpressionFile(expressionFile)

#print("X\tY")
genes=geneToFC.keys()
for gene in genes:
    enhancers=geneToEnhancer.get(gene,None)
    if(enhancers is None or len(enhancers)<1): continue
    P300=[]
    for enhancerID in enhancers:
        if(enhancerToP300.get(enhancerID,None) is None): continue
        P300.append(enhancerToP300[enhancerID])
    if(len(P300)<1): continue
    #X=max(P300)
    X=sum(P300)
    Y=geneToFC[gene]
    print(X,Y,sep="\t")
    #print(len(enhancers),file=sys.stderr)



