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
import os
import ProgramName
from Interval import Interval
from Rex import Rex
rex=Rex()
from Fastb import Fastb

MIN_PEAK_LEN=100 #250
FIX_LENGTHS=True
DELTA="/home/bmajoros/GGR/delta"
#ENHANCER_COORDS=DELTA+"/p300-hg38.bed"
PREDICTIONS=DELTA+"/genomewide-loopy.txt"
MOTIF_FASTBS=DELTA+"/motifs-continuous"
ALL_MOTIFS=set()
MOTIFS=("CTCF","AP1","GR","FOX","KLF","CEBP")
for motif in MOTIFS: ALL_MOTIFS.add(motif)

def loadPredictions(filename,wantTime):
    byChr={}
    IN=open(filename,"rt")
    for line in IN:
        fields=line.rstrip().split("\t")
        (peakID,begin,end,task,LLR,P,parse,topBottom)=fields
        if(not rex.find("(\S+)\.(t\d+)\.fastb",peakID)):
            raise Exception("can't parse"+peakID)
        peakID=rex[1]; thisTime=rex[2]
        if(thisTime!=wantTime): continue
        fields=parse.split("|")
        parse=[]
        nextElemID=1
        for field in fields:
            if(not rex.find("(\d+):(\d+)-(\d+)",field)):
                raise Exception("can't parse "+field)
            state=int(rex[1]); begin=int(rex[2]); end=int(rex[3])
            if(state!=2 and state!=3 and state!=4): continue
            if(end-begin<MIN_PEAK_LEN): continue
            if(FIX_LENGTHS):
                mid=int((begin+end)/2)
                begin=mid-int(MIN_PEAK_LEN/2)
                end=begin+MIN_PEAK_LEN
            peak=Interval(begin,end)
            peak.type="peak" if state==3 else "hump"
            peak.ID=peakID+"_"+str(nextElemID)
            nextElemID+=1
            peak.motifs={}; ### set()
            parse.append(peak)
        byChr[peakID]=parse
    IN.close()
    return byChr

def hasPeak(elements):
    for elem in elements:
        if(elem.type=="peak"): return True
    return False

def assignHitsToPeaks(motifID,hits,peaks):
    for hit in hits:
        for peak in peaks:
            if(hit.overlaps(peak)):
                #peak.motifs.add(motifID)
                overlapLen=peak.intersect(hit).getLength()
                addMotif(motifID,overlapLen,peak)

def addMotif(motifName,length,peakHump):
    peakHump.motifs[motifName]=peakHump.motifs.get(motifName,0)+length

def processFastb(filename,wantTime,parses):
    if(not rex.find("(\S+)\.standardized_across_all_timepoints\.(t\d+)\.fastb",filename)):
        raise Exception(filename)
    peakID=rex[1]; thisTime=rex[2]
    #if(thisTime!=wantTime): return ###
    peaksAndHumps=parses.get(peakID,None)
    if(peaksAndHumps is None): return
    if(not hasPeak(peaksAndHumps)): return
    fastb=Fastb(MOTIF_FASTBS+"/"+filename)
    numTracks=fastb.numTracks()
    for i in range(numTracks):
        track=fastb.getIthTrack(i)
        motifID=track.getID()
        if(motifID=="GR/AR/MR"): motifID="GR"
        hits=track.getNonzeroRegions()
        assignHitsToPeaks(motifID,hits,peaksAndHumps)
    for elem in peaksAndHumps:
        print(elem.ID,elem.type,elem.getLength(),sep="\t",end="")
        for motif in ALL_MOTIFS: #elem.motifs.keys():
            L=elem.motifs.get(motif,0)
            print("\t"+motif+"="+str(L),end="")
        print()
        elem.motifs=None

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <timepoint>\n")
(timepoint,)=sys.argv[1:]

# Load parses from RESCUE
byChr=loadPredictions(PREDICTIONS,timepoint)

# Process motif FASTB files
files=os.listdir(MOTIF_FASTBS)
for file in files:
    processFastb(file,timepoint,byChr)


