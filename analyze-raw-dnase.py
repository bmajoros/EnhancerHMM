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

CHIPS=("DNase",)
TRACK_NAME="DNase"
MOTIF_ONLY=False
MARGIN=100
CHIP_DIR="/data/reddylab/projects/GGR/data/hmm/preprocessing/EP300_not_standardized/split_fastb"
MOTIFS_DIR="/home/bmajoros/GGR/delta/motifs-continuous"

def getMotifs(trackName,peakID):
    filename=MOTIFS_DIR+"/"\
        +peakID+".standardized_across_all_timepoints.t00.fastb"
    fastb=Fastb(filename)
    track=fastb.getTrackByName(trackName)
    hits=track.getNonzeroRegions()
    return hits

def annotateMotifs(peaks,peakID):
    motifs=getMotifs(TRACK_NAME,peakID)
    for peak in peaks:
        peak.motif=False
        interval=peak
        interval.begin=interval.intCenter()-100
        interval.end=interval.intCenter()+100
        for motif in motifs:
            if(motif.overlaps(interval)): peak.motif=True

def process(peakID,parse,scores):
    peaks=getPeaks(parse)
    if(MOTIF_ONLY): annotateMotifs(peaks,peakID)
    processCHIP(peaks,peakID)
    numPeaks=len(peaks)
    for peak in peaks:
        if(MOTIF_ONLY and peak.motif==False): continue
        print(peak.score,numPeaks,peak.length(),sep="\t")
        if(scores.get(peakID,None) is None): scores[peakID]=[]
        scores[peakID].append(peak.score)

def processCHIP(peaks,peakID):
    for chip in CHIPS:
        chipfile=CHIP_DIR+"/"+chip+"."+TIME+"."+peakID+".fastb"
        fastb=Fastb(chipfile)
        getCHIPscores(peaks,fastb,chip)

def getCHIPscores(peaks,fastb,CHIP_NAME):
    track=fastb.getTrackByName(CHIP_NAME)
    for peak in peaks:
        array=track.data[peak.begin:peak.end]
        peak.score+=max(array)

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
        peak=Interval(begin,end)
        peak.score=0.0
        peaks.append(peak)
    return peaks

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
    if(rex.find("(\S+)\.(t\d+)",fastb) and rex[2]!=TIME): continue
    process(rex[1],parse,scoreSets)
IN.close()

filename="dnase-sums-motif.txt" if MOTIF_ONLY else "dnase-sums.txt"
OUT=open(filename,"wt")
for key in scoreSets.keys():
    total=sum(scoreSets[key])
    N=len(scoreSets[key])
    mean=float(total)/float(N)
    print(total,N,mean,sep="\t",file=OUT)
OUT.close()

