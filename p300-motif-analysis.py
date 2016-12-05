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
import sys
from Fastb import Fastb
from Interval import Interval

MAX_TRAJECTORIES=3000
BASE="/home/bmajoros/GGR/p300/"
MOTIF_DIR=BASE+"motif-tracks"
PARSED_DIR=BASE+"parsed"
TRAJECTORIES=BASE+"trajectory-class-pos5-neg5.txt"
TIMEPOINTS=("t00","t05","t1","t2","t3","t4","t5","t6","t7","t8","t10","t12")

peakHits={}
humpHits={}
flankHits={}
totalPeakBases=0
totalHumpBases=0
totalFlankBases=0

#=========================== FUNCTIONS =========================
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
    global totalFlankBases
    global totalPeakBases
    global totalHumpBases
    filename=PARSED_DIR+"/"+peak+".standardized_across_all_timepoints."+time+".fastb"
    fastb=Fastb(filename)
    track=fastb.getTrackByName("state")
    regions=track.getContiguousRegions()
    for region in regions:
        state=int(region.value)
        if(state==1): totalFlankBases+=analyze(region,motifs,flankHits)
        elif(state==2): totalHumpBases+=analyze(region,motifs,humpHits)
        elif(state==3): totalPeakBases+=analyze(region,motifs,peakHits)

def analyze(region,motifs,hitsHash):
    for motifName in list(motifs.keys()):
        hits=motifs[motifName]
        for hit in hits:
            if(hit.overlaps(region)):
                if(not hitsHash.get(motifName,None)): hitsHash[motifName]=0
                hitsHash[motifName]+=1
    return region.getLength()

def getRatios(hash,length):
    ratios={}
    keys=list(hash.keys())
    for key in keys:
        numHits=hash[key]
        ratio=numHits/length
        ratios[key]=ratio
    return ratios

#=========================== main() =========================
numProcessed=0
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
        numProcessed+=1
        print(numProcessed,"trajectories processed",file=sys.stderr)
        if(numProcessed>=MAX_TRAJECTORIES): break
peakHash=getRatios(peakHits,totalPeakBases)
humpHash=getRatios(humpHits,totalHumpBases)
flankHash=getRatios(flankHits,totalFlankBases)
keys=list(peakHash.keys())
print("factor\tpeak/hump\tpeak/flank\thump/flank")
for key in keys:
    # Analyze enrichment in peaks versus humps:
    peakRatio=peakHash.get(key,None)
    humpRatio=humpHash.get(key,None)
    if(peakRatio is None or humpRatio is None): exit("zero")
    normalized=int(peakRatio/humpRatio*100.0+5.0/9.0)/100.0
    # Analyze enrichment in peaks versus flanks:
    flankRatio=flankHash.get(key,None)
    if(flankRatio is None): exit("zero2")
    normalized2=int(peakRatio/flankRatio*100.0+5.0/9.0)/100.0
    normalized3=int(humpRatio/flankRatio*100.0+5.0/9.0)/100.0
    # Analyze enrichment in humps versus flanks:
    print(key+"\t"+str(normalized)+"\t"+str(normalized2)+"\t"+str(normalized3))


