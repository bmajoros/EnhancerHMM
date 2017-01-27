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
from scipy.stats import mannwhitneyu
import numpy as numpy
from Rex import Rex
rex=Rex()

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <trajectories.txt> <peaks-motifs.txt> <raw-P300.txt> <cluster-membership.txt> <which-cluster>\n")
(trajectoryFile,motifFile,p300file,clusterFile,whichCluster)=sys.argv[1:]
whichCluster=int(whichCluster)

cluster1=set()
with open(clusterFile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=2): continue
        (peak,cluster)=fields
        cluster=int(cluster)
        if(not rex.find("\"(\S+)\"",peak)): exit(peak)
        peak=rex[1]
        if(cluster==whichCluster or whichCluster==0): cluster1.add(peak)

keepPeak=set()
with open(motifFile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=2): continue
        (peak,motif)=fields
        if(motif=="GR"): keepPeak.add(peak)

numPeaks={}
with open(trajectoryFile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=12): continue
        maxPeaks=max(fields[1:])
        if(maxPeaks==0): continue
        peakID=fields[0]
        numPeaks[peakID]=int(maxPeaks)

valuesByPeaks=[]
for i in range(3): valuesByPeaks.append([])
with open(p300file,"rt") as IN:
    header=IN.readline()
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=12): continue
        peakID=fields[0]
        if(not peakID in keepPeak): continue
        if(not peakID in cluster1): continue
        values=fields[1:]
        for i in range(len(values)): values[i]=float(values[i])
        maxValue=max(values)
        sumValue=sum(values)/float(len(values))
        np=numPeaks.get(peakID,0)
        if(np==0): continue
        #valuesByPeaks[np].append(maxValue)
        valuesByPeaks[np].append(sumValue)

(U,P)=mannwhitneyu(valuesByPeaks[1],valuesByPeaks[2])
print(2**numpy.median(valuesByPeaks[1]),2**numpy.median(valuesByPeaks[2]),U,P)



