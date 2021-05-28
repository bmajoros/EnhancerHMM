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
from Interval import Interval
from SummaryStats import SummaryStats
from Rex import Rex
rex=Rex()

COORDS_FILE="p300-hg38.bed"
DISTANCES="peak-distances.txt"

def loadCoords():
    byID={}
    byChr={}
    with open(COORDS_FILE,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            (chr,begin,end,id)=fields
            if(rex.find("random",chr)): continue
            begin=int(begin); end=int(end)
            rec=Interval(begin-1000,begin+1000)
            rec.id=id
            if(byChr.get(chr,None) is None): byChr[chr]=[]
            byChr[chr].append(rec)
    for chrom in byChr.keys():
        array=byChr[chrom]
        array.sort(key=lambda x: x.begin)
        deduplicate(array)
        for rec in array: byID[rec.id]=rec
    return (byID,byChr)

def deduplicate(array):
    L=len(array)
    i=0
    while(i<L-1):
        if(array[i].overlaps(array[i+1])):
            del array[i+1]
            L-=1
        else: i+=1

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <genomewide-loopy.txt> <timepoint> <min-length>\n")
(infile,TIMEPOINT,MIN_LEN)=sys.argv[1:]
MIN_LEN=int(MIN_LEN)

(byID,byChr)=loadCoords()
singletons=0; multipeaks=0
IN=open(infile,"rt")
distances=[]
OUT=open(DISTANCES,"wt")
COUNTS=open("peak-counts.txt","wt")
for line in IN:
    fields=line.rstrip().split()
    (peakID,begin,end,taskID,LLR,posterior,parse,top)=fields
    if(not rex.find("^([^\.]+)\.(t\d+)",peakID)): raise Exception(peakID)
    peakID=rex[1]; time=rex[2]
    if(time!=TIMEPOINT): continue
    coords=byID.get(peakID,None)
    if(coords is None): continue
    subfields=parse.split("|")
    numPeaks=0
    prevPeakEnd=None
    for subfield in subfields:
        if(not rex.find("(\d+):(\d+)-(\d+)",subfield)):
            raise Exception("can't parse: "+subfield)
        state=int(rex[1]); b=int(rex[2]); e=int(rex[3])
        length=e-b
        if(length<MIN_LEN): continue
        if(state==3):
            numPeaks+=1
            if(numPeaks>1):
                dist=b-prevPeakEnd
                distances.append(dist)
                print(dist,file=OUT)
            prevPeakEnd=e
    if(numPeaks==1): singletons+=1
    elif(numPeaks>1): multipeaks+=1
    print(numPeaks,file=COUNTS)
IN.close()
OUT.close()
COUNTS.close()
percent=round(100.0*float(multipeaks)/float(singletons+multipeaks),2)
print(singletons,"singletons,",multipeaks,"multipeaks = ",percent,"%")
(mean,SD,min,max)=SummaryStats.summaryStats(distances)
print("distances: ",mean,"+/-",SD,"min=",min,"max=",max)

