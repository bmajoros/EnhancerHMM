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
from BedReader import BedReader
from Bed3Record import Bed3Record
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
        #print(distance,nearest.interval.center(),elem.interval.center(),nearest.chr,elem.chr,sep="\t")
        print(elem.numPeaks,distance,sep="\t")

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

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <enhancer-locations.bed> <DEGs.bed>\n")
(enhancerFile,degFile)=sys.argv[1:]

enhancerHash=BedReader.hashBySubstrate(enhancerFile)
degHash=BedReader.hashBySubstrate(degFile)
chroms=enhancerHash.keys()
for chr in chroms:
    if(len(chr)>5): continue
    enhancers=enhancerHash[chr]
    labelEnhancers(enhancers)
    degs=degHash[chr]
    labelGenes(degs)
    combined=enhancers
    combined.extend(degs)
    combined.sort(key=lambda x: x.interval.center())
    getDistances(combined)



