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

MARGIN=100
P300_DIR="/data/reddylab/projects/GGR/subprojects/hmm/data/EP300_not_standardized"

def process(peakID,parse,scores):
    p300file=P300_DIR+"/"+peakID+".fastb"
    peaks=getPeaks(parse)
    fastb=Fastb(p300file)
    getP300scores(peaks,fastb,peakID)
    numPeaks=len(peaks)
    for peak in peaks:
        print(peak.score,numPeaks,peak.length(),sep="\t")
        if(scores.get(peakID,None) is None): scores[peakID]=[]
        scores[peakID].append(peak.score)

def getPeaks(parse):
    peaks=[]
    fields=parse.split("|")
    for field in fields:
        if(not rex.find("(\d+):(\d+)-(\d+)",field)):
            exit("can't parse field: "+field)
        state=int(rex[1])
        if(state!=3): continue
        begin=int(rex[2]); end=int(rex[3])
        peakLen=end-begin
        if(peakLen<MIN_PEAK_LEN): continue
        if(MARGIN>0):
            mid=int((begin+end)/2)
            begin=mid-MARGIN
            end=mid+MARGIN
        if(begin<0 or end>2000): continue
        peaks.append(Interval(begin,end))
    return peaks

def getP300scores(peaks,fastb,peakID):
    ###
    if(not rex.find("(\S+)\.t\d+",peakID)): exit(peakID)
    id=rex[1]
    otherFastb=Fastb(P300_DIR+"/"+id+".t00.fastb")
    otherTrack=otherFastb.getTrackByName("EP300")
    ###
    track=fastb.getTrackByName("EP300")
    for peak in peaks:
        array=track.data[peak.begin:peak.end]
        otherArray=otherTrack.data[peak.begin:peak.end]
        #peak.score=max(array)
        peak.score=max(array)-max(otherArray)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()
         +" <genomewide.txt> <LLR-threshold> <min-peak-len> <time>\n")
(infile,threshold,MIN_PEAK_LEN,TIME)=sys.argv[1:]
MIN_PEAK_LEN=int(MIN_PEAK_LEN)

scoreSets={}
IN=open(infile,"rt")
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=8): continue
    (fastb,b,e,taskID,LLR,P,parse,top)=fields
    if(LLR<threshold): continue
    if(rex.find("(\S+)\.fastb",fastb)): fastb=rex[1]
    if(rex.find("\S+\.(t\d+)",fastb) and rex[1]!=TIME): continue
    process(fastb,parse,scoreSets)
IN.close()

OUT=open("p300-sums.txt","wt")
for key in scoreSets.keys():
    total=sum(scoreSets[key])
    N=len(scoreSets[key])
    mean=float(total)/float(N)
    print(total,N,mean,sep="\t",file=OUT)
OUT.close()

