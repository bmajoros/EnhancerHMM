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
from scipy.stats import fisher_exact

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <motifs-peaks-humps.txt>\n")
(motifsFile,)=sys.argv[1:]

numPeaks=0; numHumps=0
peakLength=0; humpLength=0
peakCounts={}; humpCounts={}
with open(motifsFile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        type=fields[1]
        if(type=="peak"):
            numPeaks+=1
            peakLength+=float(fields[2])
            for field in fields[3:]:
                peakCounts[field]=peakCounts.get(field,0)+1
        elif(type=="hump"):
            numHumps+=1
            humpLength+=float(fields[2])
            for field in fields[3:]:
                humpCounts[field]=humpCounts.get(field,0)+1
        else: raise Exception(type)
motifs=set()
for motif in peakCounts.keys(): motifs.add(motif)
for motif in humpCounts.keys(): motifs.add(motif)
array=[]
for motif in motifs:
    peakCount=peakCounts.get(motif,0)
    humpCount=humpCounts.get(motif,0)
    peakFreq=round(1000000*float(peakCount)/float(peakLength),1)
    humpFreq=round(1000000*float(humpCount)/float(humpLength),1)
    enrichment=round(peakFreq/humpFreq,1)
    #print(motif,enrichment,peakFreq,humpFreq,peakCount,humpCount,sep="\t")
    peakAbsent=numPeaks-peakCount
    humpAbsent=numHumps-humpCount
    peakExpected=int(float(humpCount)/float(numHumps)*numPeaks)
    enrichment=float(peakCount-peakExpected)/float(peakExpected)
    (oddsRatio,P)=fisher_exact([[peakCount,peakAbsent],[humpCount,humpAbsent]])
    m=motif
    if(m=="GR/AR/MR"): m="GR"
    array.append([m,P,enrichment,peakCount,peakAbsent,humpCount,humpAbsent])
array.sort(key=lambda x: -x[2])
for rec in array:
    (motif,P,enrichment,peakCount,peakAbsent,humpCount,humpAbsent)=rec
    print(motif,P,round(enrichment,2),\
              peakCount,peakAbsent,humpCount,humpAbsent,sep="\t")
