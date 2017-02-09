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
from Fastb import Fastb
from Rex import Rex
rex=Rex()

P300_DIR="/data/reddylab/projects/GGR/subprojects/hmm/data/EP300_not_standardized"

def process(peakID,parse):
    p300file=P300_DIR+"/"+peakID+".fastb"
    peaks=getPeaks(parse)
    fastb=Fastb(p300file)
    getP300scores(peaks,fastb)
    numPeaks=len(peaks)
    for peak in peaks:
        print(peak.score,numPeaks,peak.length(),sep="\t")

def getPeaks(parse):
    peaks=[]
    fields=parse.split("|")
    for field in fields:
        if(not rex.find("(\d+):(\d+)-(\d+)",field)):
            exit("can't parse field: "+field)
        state=int(rex[1])
        if(state!=3): continue
        begin=int(rex[2]); end=int(rex[3])
        peaks.append(Interval(begin,end))
    return peaks

def getP300scores(peaks,fastb):
    track=fastb.getTrackByName("EP300")
    for peak in peaks:
        array=track.data[peak.begin:peak.end]
        peak.score=max(array)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <genomewide.txt> <LLR-threshold>\n")
(infile,threshold)=sys.argv[1:]

IN=open(infile,"rt")
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=8): continue
    (fastb,b,e,taskID,LLR,P,parse,top)=fields
    if(LLR<threshold): continue
    if(rex.find("(\S+)\.fastb",fastb)): fastb=rex[1]
    if(rex.find("\S+\.(t\d+)",fastb) and rex[1]!="t05"): continue
    process(fastb,parse)
IN.close()



