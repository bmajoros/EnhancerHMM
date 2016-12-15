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
from Interval import Interval
from Rex import Rex
rex=Rex()

def readPath(inFile):
    array=[]
    with open(inFile,"rt") as IN:
        header=IN.readline()
        for line in IN:
            if(not rex.find("\S+",line)): continue
            state=int(line)
            array.append(state)    
    return array

def findElements(path):
    elements=[]
    begin=-1
    L=len(path)
    for i in range(L):
        if(i==0 and path[i]!=1 or i>0 and path[i]!=1 and path[i-1]==1):
            begin=i
        elif(i>0 and path[i-1]!=1 and path[1]==1):
            elements.append(Interval(begin,i))
            elements[len(elements)-1].states=path[begin:i]
    if(L>0 and path[L-1]!=1):
        elements.append(Interval(begin,L))
        elements[len(elements)-1].states=path[begin:L]
    return elements

def getSections(states):
    sections=[]
    L=len(states)
    elem=[]
    for i in range(L):
        state=states[i]
        eLen=len(elem)
        if(eLen==0 or elem[eLen-1]!=state): elem.append(state)
        else:
            sections.append(elem)
            elem=[]
    return sections

def satisfies(elem,minHump,minPeak,minTotal):
    sections=getSections(elem.states)
    for section in sections:
        L=len(section)
        if(L<minTotal): return False
        type=section[0]
        if(type==2 and L<minHump): return False
        if(type==3 and L<minPeak): return False
        if(type==4 and L<minHump): return False
    return True

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=7):
    exit(ProgramName.get()+
         " <chr> <in.path> <min-hump> <min-peak> <min-total> <out.bed>\n")
(chr,inFile,minHump,minPeak,minTotal,outFile)=sys.argv[1:]
minHump=int(minHump); minPeak=int(minPeak); minTotal=int(minTotal)

path=readPath(inFile)
elements=findElements(path)
for elem in elements:
    if(satisfies(elem,minHump,minPeak,minTotal)):
        print(chr+"\t"+str(elem.begin)+"\t"+str(elem.end))


