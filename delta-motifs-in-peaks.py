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
from Interval import Interval
from Fastb import Fastb
from FastbTrack import FastbTrack
from Rex import Rex
rex=Rex()

MIN_PEAK_LEN=20
MIN_SPACER_LEN=20

def getMotifs(filename):
    fastb=Fastb(filename)
    hits=[]
    for i in range(fastb.numTracks()):
        track=fastb.getIthTrack(i)
        if(not track.isDiscrete()): continue
        motifID=track.getID()
        if(motifID=="GR/AR/MR"): motifID="GR"
        getMotifIntervals(track,motifID,hits)
    return hits
    pass

def getMotifIntervals(track,id,into):
    intervals=track.getDiscreteContiguousRegions()
    for interval in intervals:
        if(interval.value=="0"): continue
        interval.motif=id
        into.append(interval)

def getMotifsInPeak(motifs,peak):
    array=[]
    for motif in motifs:
        if(motif.overlaps(peak)): array.append(motif)
    return array

def emitPairs(peakID,motifsPeak1,motifsPeak2):
    print("peak "+peakID)
    seen=set()
    for motif1 in motifsPeak1:
        for motif2 in motifsPeak2:
            id1=motif1.motif; id2=motif2.motif
            if(id1>id2): (id1,id2)=(id2,id1)
            key=id1+"\t"+id2
            if(key in seen): continue
            seen.add(key)
            print("\t"+key)
    pass

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <genomewide.txt> <motifs-discrete-dir>\n")
(peaksFile,motifsDir)=sys.argv[1:]

IN=open(peaksFile,"rt")
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=7): continue
    (fastb,begin,end,id,LLR,parse,topBottom)=fields
    if(not rex.find("^(\S+)\.t",fastb)): exit("can't parse "+fastb)
    peakID=rex[1]
    if(not rex.find("(\d+)-(\d+):(\d+)-(\d+):(\d+)-(\d+)",parse)):
        exit("can't parse path: "+parse)
    peak1=Interval(int(rex[1]),int(rex[2]))
    peak2=Interval(int(rex[5]),int(rex[6]))
    spacer=Interval(int(rex[3]),int(rex[4]))
    if(peak1.length()<MIN_PEAK_LEN or peak2.length()<MIN_PEAK_LEN or
       spacer.length()<MIN_SPACER_LEN): continue
    motifFile=motifsDir+"/"+peakID+\
        ".standardized_across_all_timepoints.t00.fastb"
    motifs=getMotifs(motifFile)
    motifsPeak1=getMotifsInPeak(motifs,peak1)
    motifsPeak2=getMotifsInPeak(motifs,peak2)
    emitPairs(peakID,motifsPeak1,motifsPeak2)
IN.close()


