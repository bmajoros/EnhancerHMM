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
from BedReader import BedReader
from Interval import Interval
from Rex import Rex
rex=Rex()

MIN_PEAK_LEN=250
STANDARD_PEAK_WIDTH=200 # 200
DELTA="/home/bmajoros/GGR/delta"
KEVIN=DELTA+"/kevin-intersection.bed"
ENHANCER_COORDS=DELTA+"/p300-hg38.bed"
PREDICTIONS=DELTA+"/genomewide-loopy.txt"

def loadEnhancerCoords(filename):
    hash={}
    list=BedReader.readAll(filename)
    for elem in list:
        hash[elem.name]=elem
    return hash

def loadPeaks(filename,enhancerCoords,whichTime):
    byChr={}
    with open(filename) as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=8): continue
            (fastb,b,e,taskID,LLR,P,parse,top)=fields
            #if(float(LLR)<1000): continue ###
            if(not rex.find("(\S+)\.(t\d+)\.fastb",fastb)):
                raise Exception(fastb)
            enhancerID=rex[1]; time=rex[2]
            if(time!=whichTime): continue
            chunkCoords=enhancerCoords.get(enhancerID,None)
            origin=chunkCoords.interval.begin
            #if(chunkCoords is None): continue
            subfields=parse.split("|")
            peaks=[]
            numPeaks=0
            for subfield in subfields:
                if(not rex.find("(\d+):(\d+)-(\d+)",subfield)):
                    raise Exception(subfield)
                state=int(rex[1])
                if(state!=2 and state!=3 and state!=4): continue
                peaktype="peak" if state==3 else "hump"
                if(peaktype=="peak"): numPeaks+=1
                begin=int(rex[2]); end=int(rex[3])
                peakLen=end-begin
                if(peakLen<MIN_PEAK_LEN): continue
                peak=Interval(origin+begin,origin+end)
                peak.type=peaktype
                peak.id=enhancerID+"_"+str(numPeaks-1)
                if(STANDARD_PEAK_WIDTH>0):
                    peak.begin=peak.intCenter()-STANDARD_PEAK_WIDTH/2
                    peak.end=peak.intCenter()+STANDARD_PEAK_WIDTH/2
                peak.chr=chunkCoords.chr
                peak.motifs=[]
                peaks.append(peak)
            #numPeaks=len(peaks)
            if(numPeaks<1): continue
            #enhancerType="singleton" if numPeaks==1 else "multipeak"
            #for peak in peaks: peak.type=enhancerType
            if(byChr.get(chunkCoords.chr,None) is None):
                byChr[chunkCoords.chr]=[]
            for peak in peaks: byChr[chunkCoords.chr].append(peak)
    return byChr

def loadTFs(filename,hash):
    list=BedReader.readAll(filename)
    for tf in list:
        array=hash.get(tf.chr,None)
        if(array is None): continue
        rec=tf.interval
        rec.name=tf.name
        rec.type="TF"
        array.append(rec)

def sortChroms(hash):
    chroms=hash.keys()
    for chr in chroms:
        array=hash[chr]
        array.sort(key=lambda x: x.begin)

def putTFsInPeaks(peaksByChr):
    chroms=peaksByChr.keys()
    for chr in chroms:
        array=peaksByChr[chr]
        L=len(array)
        for i in range(L):
            peak=array[i]
            if(peak.type!="peak" and peak.type!="hump"): continue
            #if(peak.type!="singleton" and peak.type!="multipeak"): continue
            scanLeft(i,array,peak)
            scanRight(i,array,L,peak)
        newArray=[]
        for elem in array:
            if(elem.type!="peak" and elem.type!="hump"): continue
            #if(elem.type=="singleton" or elem.type=="multipeak"):
            newArray.append(elem)
        peaksByChr[chr]=newArray

def scanLeft(index,array,peak):
    for i in range(index-1,-1,-1):
        tf=array[i]
        if(tf.type!="TF"): continue
        if(tf.overlaps(peak)): peak.motifs.append(tf)
        elif(tf.distance(peak)>30): break

def scanRight(index,array,L,peak):
    for i in range(index+1,L):
        tf=array[i]
        if(tf.type!="TF"): continue
        if(tf.overlaps(peak)): peak.motifs.append(tf)
        elif(tf.distance(peak)>30): break

def dump(peaksByChr):
    chroms=peaksByChr.keys()
    for chr in chroms:
        array=peaksByChr[chr]
        for peak in array:
            print(peak.id,peak.type,peak.getLength(),sep="\t",end="")
            for tf in peak.motifs:
                print("\t"+tf.name,end="")
            print()

#=========================================================================
# main()
#=========================================================================

if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <timepoint>")
(TIMEPOINT,)=sys.argv[1:]

enhancerIdToCoords=loadEnhancerCoords(ENHANCER_COORDS)
peaksByChr=loadPeaks(PREDICTIONS,enhancerIdToCoords,TIMEPOINT)
loadTFs(KEVIN,peaksByChr)
sortChroms(peaksByChr)
putTFsInPeaks(peaksByChr)
dump(peaksByChr)



