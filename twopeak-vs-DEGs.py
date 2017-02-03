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
import random
import ProgramName
from BedReader import BedReader
from Bed3Record import Bed3Record
from SummaryStats import SummaryStats
from Rex import Rex
rex=Rex()

def labelEnhancers(elements):
    for elem in elements:
        if(not rex.find("(\S+)_(\d+)",elem.name)):
            exit("can't parse "+elem.name)
        elem.tag="enhancer"
        elem.name=rex[1]
        elem.numPeaks=2 if "2" in rex[2] else 1

def labelGenes(elements):
    for elem in elements:
        elem.tag="gene"

def getDistances(elements):
    for i in range(len(elements)):
        elem=elements[i]
        if(elem.tag!="enhancer"): continue
        nearest=findNearestGene(elements,i)
        if(nearest is None): continue
        distance=int(abs(nearest.interval.center()-elem.interval.center()))
        print(elem.numPeaks,distance,sep="\t")

def getDistancesSingletons(elements):
    distances=[]
    for i in range(len(elements)):
        elem=elements[i]
        if(elem.tag!="enhancer"): continue
        if(elem.numPeaks!=1): continue
        nearest=findNearestGene(elements,i)
        if(nearest is None): continue
        distance=int(abs(nearest.interval.center()-elem.interval.center()))
        distances.append(distance)
    return distances

def findNearestGene(elements,to):
    left=searchLeft(elements,to)
    right=searchRight(elements,to)
    if(left is None and right is None): return None
    if(left is None): return right
    if(right is None): return left
    leftDist=elements[to].interval.center()-left.interval.center()
    rightDist=right.interval.center()-elements[to].interval.center()
    if(leftDist<rightDist): return left
    return right

def searchLeft(elements,to):
    for i in range(to-1,-1,-1):
        if(elements[i].tag=="gene"): return elements[i]
    return None

def searchRight(elements,to):
    for i in range(to+1,len(elements)):
        if(elements[i].tag=="gene"): return elements[i]
    return None

def randomize(array):
    N=len(array)
    L=array[N-1].interval.end
    for elem in array:
        if(elem.tag=="enhancer"):
            elem.interval.begin=int(random.randint(1,L))
            elem.interval.end=elem.interval.begin+1
    combined.sort(key=lambda x: x.interval.center())

def thin(array1,array2):
    L1=len(array1)
    L2=len(array2)
    if(L1>L2): thinArray(array1,L2)
    elif(L2>L1): thinArray(array2,L1)

def shuffle(array):
    L=len(array)
    for i in range(L):
        j=random.randint(0,L-1)
        a=array[i]
        b=array[j]
        array[i]=b
        array[j]=a

def thinArray(array,size):
    shuffle(array)
    del array[size:]

def getByPeakCount(array,count):
    ret=[]
    for elem in array:
        if(elem.numPeaks==count): ret.append(elem)
    return ret

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <enhancer-locations.bed> <DEGs.bed> <randomize:0|1>\n")
(enhancerFile,degFile,wantRandomize)=sys.argv[1:]
wantRandomize=int(wantRandomize)!=0

NUM_ITER=1000

for i in range(NUM_ITER):
    enhancerHash=BedReader.hashBySubstrate(enhancerFile)
    degHash=BedReader.hashBySubstrate(degFile)
    chroms=enhancerHash.keys()
    for chr in chroms:
        if(len(chr)>5): continue
        enhancers=enhancerHash[chr]
        labelEnhancers(enhancers)

        # Thinning:
        singlePeaks=getByPeakCount(enhancers,1)
        dualPeaks=getByPeakCount(enhancers,2)
        thin(singlePeaks,dualPeaks)
        enhancers=singlePeaks
        enhancers.extend(dualPeaks)

        degs=degHash[chr]
        labelGenes(degs)
        combined=enhancers
        combined.extend(degs)
        combined.sort(key=lambda x: x.interval.center())
        if(wantRandomize): randomize(combined)
        #getDistances(combined)
        distances=getDistancesSingletons(combined)
        [mean,SD,min,max]=SummaryStats.summaryStats(distances)
        print(mean)



