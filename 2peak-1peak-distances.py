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
import math
import ProgramName
from BedReader import BedReader
from Bed3Record import Bed3Record
from Rex import Rex
rex=Rex()

FROM_PEAKS=2
TO_PEAKS=1

def labelElements(elements):
    for elem in elements:
        if(not rex.find("(\S+)_(\d+)",elem.name)):
            exit("can't parse "+elem.name)
        elem.name=rex[1]
        elem.numPeaks=2 if rex[2]=="20" else 1

def getDistances(elements):
    for i in range(len(elements)):
        elem=elements[i]
        if(elem.numPeaks!=FROM_PEAKS): continue
        nearest=findNearest(elements,i,TO_PEAKS)
        if(nearest is None): continue
        distance=int(abs(nearest.interval.center()-elem.interval.center()))
        #print(distance,nearest.interval.center(),elem.interval.center(),nearest.chr,elem.chr,sep="\t")
        print(distance)

def findNearest(elements,to,numPeaks):
    left=searchLeft(elements,to,numPeaks)
    right=searchRight(elements,to,numPeaks)
    #print("left",left,"right",right)
    if(left is None and right is None): return None
    if(left is None): return right
    if(right is None): return left
    leftDist=elements[to].interval.center()-left.interval.center()
    rightDist=right.interval.center()-elements[to].interval.center()
    if(leftDist<rightDist): return left
    return right

def searchLeft(elements,to,numPeaks):
    for i in range(to-1,-1,-1):
        if(elements[i].numPeaks==numPeaks): return elements[i]
    return None

def searchRight(elements,to,numPeaks):
    for i in range(to+1,len(elements)):
        if(elements[i].numPeaks==numPeaks): return elements[i]
    return None

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <peak-locations.bed>\n")
infile=sys.argv[1]

bySubstrate=BedReader.hashBySubstrate(infile)
chroms=bySubstrate.keys()
for chr in chroms:
    elements=bySubstrate[chr]
    elements.sort(key=lambda x: x.interval.center())
    labelElements(elements)
    getDistances(elements)


