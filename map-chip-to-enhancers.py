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

NEXT_ID=1

def findLeft(array,index,enhancer,found):
    for i in range(index,-1,-1):
        if(array[i].interval.overlaps(enhancer.interval) and
           array[i].type=="chip"):
            found.append(array[i])
        elif(enhancer.interval.begin-array[i].interval.end>1000): break
    return found

def findRight(array,index,enhancer,found):
    for i in range(index,len(array)):
        if(array[i].interval.overlaps(enhancer.interval) and
           array[i].type=="chip"):
            found.append(array[i])
        elif(array[i].interval.begin-enhancer.interval.end>1000): break
    return found

def getOverlappingChips(array,index,enhancer):
    found=[]
    findLeft(array,index,enhancer,found)
    findRight(array,index,enhancer,found)
    return found

def mapChip(array):
    global NEXT_ID
    L=len(array)
    for i in range(L):
        enhancer=array[i]
        if(enhancer.type!="enhancer"): continue
        chips=getOverlappingChips(array,i,enhancer)
        for chip in chips:
            local=chip.interval.intersect(enhancer.interval)
            local.shift(-enhancer.interval.begin)
            name="elem"+str(NEXT_ID)
            NEXT_ID+=1
            score=0
            print(enhancer.name,local.begin,local.end,name,score,sep="\t")

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <chip.bed> <p300-hg19.bed>\n")
(chipBed,enhancerBed)=sys.argv[1:]

chipHash=BedReader.hashBySubstrate(chipBed)
enhancerHash=BedReader.hashBySubstrate(enhancerBed)

for chr in enhancerHash.keys():
    if(len(chr)>5): continue
    array=chipHash[chr]
    for elem in array: elem.type="chip"
    array2=enhancerHash[chr]
    for elem in array2:
        elem.type="enhancer"
        elem.interval.begin=elem.interval.intCenter()-1000
        elem.interval.end=elem.interval.intCenter()+1000
    array.extend(array2)
    array.sort(key=lambda x: x.interval.begin)
    mapChip(array)
