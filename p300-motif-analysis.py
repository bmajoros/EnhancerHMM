#!/bin/env python
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
from Fastb import Fastb
from Interval import Interval

BASE="/home/bmajoros/GGR/p300/"
MOTIF_DIR=BASE+"motif-tracks"
PARSED_DIR=BASE+"parsed"
TRAJECTORIES=BASE+"trajectory-class-pos5-neg5.txt"
TIMEPOINTS=("t00","t05","t1","t2","t3","t4","t5","t6","t7","t8","t10","t12")

peakHits={}
humpHits={}
totalPeakBases=0
totalHumpBases=0

def loadMotifs(peak):
    filename=MOTIF_DIR+"/"+peak+".standardized_across_all_timepoints.t00.fastb"
    fastb=Fastb(filename)
    motifs={}
    for i in range(fastb.numTracks()):
        track=fastb.getIthTrack(i)
        hits=track.getNonzeroRegions()
        if(len(hits)>0): motifs[track.getID()]=hits
    return motifs
        

def parse(peak,time,motifs):
    filename=PARSED_DIR+"/"+peak+".standardized_across_all_timepoints."+time+".fastb"
    fastb=Fastb(filename)
    track=fastb.getTrackByName("state")
    regions=track.getContiguousRegions()
    for region in regions:
        state=int(region.value)
        if(state==2): analyzeHump(region,motifs)
        elif(state==3): analyzePeak(region.motifs)

def analyzeHump(region,motifs):
    totalHumpBases+=region.getLength()
    for motifName in list(motifs.keys()):
        hits=motifs[motifName]
        for hit in hits:
            if(hit.overlaps(region)):
                humpHits[motifName]+=1

def analyzePeak(region,motifs):
    totalPeakBases+=region.getLength()
    for motifName in list(motifs.keys()):
        hits=motifs[motifName]
        for hit in hits:
            if(hit.overlaps(region)):
                peakHits[motifName]+=1

def getRatios(hash,length,label):
    keys=list(hash.keys())
    for key in keys:
        numHits=hash[key]
        ratio=numHits/length
        print(label,numHits,length,ratio)

with open(TRAJECTORIES,"rt") as TRAJ:
    for line in TRAJ:
        fields=line.split()
        if(len(fields)<14): continue
        peak=fields.pop(0)
        shape=fields.pop(0)
        times=[]
        for i in range(12):
            label=fields[i]
            if(label!="1"): continue
            time=TIMEPOINTS[i]
            times.append(time)
        if(len(times)==0): continue
        motifs=loadMotifs(peak)
        for time in times:
            parse(peak,time,motifs)
getRatios(peakHits,totalPeakBases,"peaks")
getRatios(humpHits,totalHumpBases,"humps")
