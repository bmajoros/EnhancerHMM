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
from Rex import Rex
rex=Rex()

TIMEPOINTS=("t05","t1","t2","t3","t4","t5","t6","t7","t8","t10","t12")

def getTimeIndex(time):
    for i in range(len(TIMEPOINTS)):
        if(time==TIMEPOINTS[i]): return i
    raise Exception("time "+time+" not found")

if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <genomewide.txt> <threshold>\n")
(infile,THRESHOLD)=sys.argv[1:]
THRESHOLD=float(THRESHOLD)

peaks={}
IN=open(infile,"rt")
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=7): continue
    (fastb,begin,end,id,LLR,parse,topBottom)=fields
    LLR=float(LLR)
    numPeaks=0
    if(LLR>=THRESHOLD):
        numPeaks=1 if topBottom=="bottom" else 2
    if(not rex.find("(\S+)\.(t\d+)\.fastb",fastb)): exit("can't parse "+fastb)
    peak=rex[1]; time=rex[2]
    timeIndex=getTimeIndex(time)
    if(peaks.get(peak,None) is None): peaks[peak]=[0]*11
    peaks[peak][timeIndex]=numPeaks
IN.close()
for peak in peaks.keys():
    trajectory=peaks[peak]
    print(peak+"\t",end="")
    for i in range(len(trajectory)):
        print(trajectory[i],end="\t" if i+1<len(trajectory) else "\n")


