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

MARGIN=250 # half width
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
        if(ONLY_MOTIF!="ALL" and motifID!=ONLY_MOTIF): continue
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

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <genomewide.txt> <motifs-discrete-dir> <motif-or-ALL> <1peak.out> <2peak.out>\n")
(peaksFile,motifsDir,ONLY_MOTIF,outfile1,outfile2)=sys.argv[1:]

OUT1=open(outfile1,"wt")
OUT2=open(outfile2,"wt")
IN=open(peaksFile,"rt")
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=7): continue
    (fastb,begin,end,id,LLR,parse,topBottom)=fields
    if(float(LLR)<1000): continue
    if(not rex.find(".t05",fastb)): continue
    if(not rex.find("^(\S+)\.t",fastb)): exit("can't parse "+fastb)
    peakID=rex[1]
    if(not rex.find("(\d+)-(\d+):(\d+)-(\d+):(\d+)-(\d+)",parse)):
        exit("can't parse path: "+parse)
    motifFile=motifsDir+"/"+peakID+\
        ".standardized_across_all_timepoints.t00.fastb"
    motifs=getMotifs(motifFile)
    numPeaks=0
    if(topBottom=="top"): numPeaks=2
    elif(topBottom=="bottom"): numPeaks=1
    if(numPeaks==2):
        peak1=Interval(int(rex[1]),int(rex[2]))
        spacer=Interval(int(rex[3]),int(rex[4]))
        peak2=Interval(int(rex[5]),int(rex[6]))
        peak1Len=peak1.length(); peak2Len=peak2.length()
        if(spacer.length()<MIN_SPACER_LEN or peak1Len<MIN_PEAK_LEN or
           peak2Len<MIN_PEAK_LEN): continue
        peak1=Interval(peak1.intCenter()-MARGIN,peak1.intCenter()+MARGIN)
        peak2=Interval(peak2.intCenter()-MARGIN,peak2.intCenter()+MARGIN)
        motifsPeak1=getMotifsInPeak(motifs,peak1)
        motifsPeak2=getMotifsInPeak(motifs,peak2)
        if(peak1.begin<0 or peak1.end>=2000 or 
           peak2.begin<0 or peak2.end>=2000): continue
        print(len(motifsPeak1),peak1Len,peak1.begin,peak1.end,
              sep="\t",file=OUT2)
        print(len(motifsPeak2),peak2Len,peak2.begin,peak2.end,
              sep="\t",file=OUT2)
    elif(numPeaks==1):
        peak=Interval(int(rex[3]),int(rex[4]))
        peakLen=peak.length()
        if(peakLen<MIN_PEAK_LEN): continue
        peak=Interval(peak.intCenter()-MARGIN,peak.intCenter()+MARGIN)
        if(peak.begin<0 or peak.end>=2000): continue
        motifs=getMotifsInPeak(motifs,peak)
        print(len(motifs),peakLen,peak.begin,peak.end,sep="\t",file=OUT1)
IN.close(); OUT1.close(); OUT2.close()


