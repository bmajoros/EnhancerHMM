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
from SummaryStats import SummaryStats
from Rex import Rex
rex=Rex()

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <trajectories.txt> <enhancer-strength.txt> <out-1peak.txt> <out-2peaks.txt>\n")
(trajectoryFile,strengthFile,outfile1,outfile2)=sys.argv[1:]

trajectories={}
with open(trajectoryFile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=12): continue
        id=fields[0]
        n=int(max(fields[1:]))
        trajectories[id]=n

OUT1=open(outfile1,"wt")
OUT2=open(outfile2,"wt")
logFC={}
FC={}
with open(strengthFile,"rt") as IN:
    header=IN.readline()
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=3): continue
        (id,log_fc,fdr)=fields
        log_fc=float(log_fc)
        if(trajectories.get(id,None) is None): continue
        numPeaks=trajectories[id]
        if(numPeaks==0): continue
        if(logFC.get(numPeaks,None) is None): logFC[numPeaks]=[]
        logFC[numPeaks].append(log_fc)
        fc=2**log_fc
        if(FC.get(numPeaks,None) is None): FC[numPeaks]=[]
        FC[numPeaks].append(fc)
        print(log_fc,file=OUT1 if numPeaks==1 else OUT2)
OUT1.close(); OUT2.close()

print("#peaks\tmeanLogFC\tSD\tmeanFC\tSD\tN")
for numPeaks in range(1,3):
    (meanFC,SDfc,minFC,maxFC)=SummaryStats.roundedSummaryStats(FC[numPeaks])
    (meanLogFC,SDlogFC,minLogFC,maxLogFC)=\
        SummaryStats.roundedSummaryStats(logFC[numPeaks])
    N=len(FC[numPeaks])
    print(numPeaks,meanLogFC,SDlogFC,meanFC,SDfc,N,sep="\t")



