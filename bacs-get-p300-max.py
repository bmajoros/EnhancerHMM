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
import os
import sys
import math
from BedReader import BedReader
from Pipe import Pipe
import SumLogProbs
from Interval import Interval
import ProgramName
import TempFilename
from Fastb import Fastb
from Rex import Rex
rex=Rex()

if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <t00 or t3> <min-elem-length>\n")
(timepoint,minLen)=sys.argv[1:]
minLen=int(minLen)

MIN_LLR=100
minHump=1
minPeak=1
fgPrior=math.log(0.5)
bgPrior=math.log(0.5)
BASE="/home/bmajoros/GGR/bacs"
P300_BED=BASE+"/p300-bacs-white.bed"
MUMMIE=os.environ["MUMMIE"]

# P300 only:
fgHMM=BASE+"/hmm/p300-only.hmm"
bgHMM=BASE+"/hmm/bg-p300-only.hmm"
fastbDir=BASE+"/fastb-p300-only"
SCHEMA=BASE+"/hmm/p300-only.schema"

# Full model:
#fgHMM=BASE+"/hmm/trained-pos5.hmm"
#bgHMM=BASE+"/hmm/trained-neg1.hmm"
#fastbDir=BASE+"/fastb"
#SCHEMA=BASE+"/hmm/p300.schema"

def sliceFastb(record,timepoint,outfile):
    begin=int(record.interval.begin)
    end=int(record.interval.end)
    #if(end-begin<400):
    #    center=int((begin+end)/2)
    #    begin=center-HALF_WIDTH
    #    end=center+HALF_WIDTH
    cmd=MUMMIE+"/fastb-slice "+SCHEMA+" "+fastbDir+"/"+record.chr+"."+timepoint+".fastb "+str(begin)+" "+str(end)+" "+outfile
    os.system(cmd)

def getLL(fastb,hmm):
    cmd=MUMMIE+"/get-likelihood "+hmm+" "+fastb
    pipe=Pipe(cmd)
    line=pipe.readline()
    return float(line)

def getPosterior(fastb):
    fgLL=getLL(fastb,fgHMM)
    bgLL=getLL(fastb,bgHMM)
    jointBG=bgLL+bgPrior
    jointFG=fgLL+fgPrior
    denom=SumLogProbs.sumLogProbs2(jointBG,jointFG)
    posterior=math.exp(jointFG-denom)
    return posterior

def getLLR(fastb):
    fgLL=getLL(fastb,fgHMM)
    bgLL=getLL(fastb,bgHMM)
    return fgLL-bgLL

def getPath(hmm,fastb):
    array=[]
    cmd=MUMMIE+"/parse "+hmm+" "+fastb
    pipe=Pipe(cmd)
    pipe.readline() # header
    while(True):
        line=pipe.readline()
        if(not line): break
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

def getSections(path):
    sections=[]
    L=len(path)
    elem=[]
    for i in range(L):
        state=path[i]
        eLen=len(elem)
        if(eLen==0 or elem[eLen-1]==state): elem.append(state)
        else:
            sections.append(elem)
            elem=[]
    return sections

def satisfiesConstraints(path,minHump,minPeak,minTotal):
    sections=getSections(path)
    totalLen=0
    for section in sections:
        L=len(section)
        type=section[0]
        if(type==2 and L<minHump): return False
        if(type==3 and L<minPeak): return False
        if(type==4 and L<minHump): return False
        if(type>1 and type<5): totalLen+=L
    if(totalLen<minTotal): return False
    return True

def applyConstraints(hmm,fastb,minHump,minPeak,minLen):
    path=getPath(hmm,fastb)
    return satisfiesConstraints(path,minHump,minPeak,minLen)

def getMax(fastbFile):
    fastb=Fastb(fastbFile)
    track=fastb.getTrackByName("EP300")
    score=track.getMax()
    return score

#=========================================================================
# main()
#=========================================================================
tempFile=TempFilename.generate(".fastb")
nextID=1
reader=BedReader(P300_BED)
while(True):
    record=reader.nextRecord()
    if(not record): break
    sliceFastb(record,timepoint,tempFile)
    id="elem"+str(nextID)
    nextID+=1
    p300=getMax(tempFile)
    print(record.chr+"\t"+str(record.interval.begin)+"\t"+str(record.interval.end)+"\t"+id+"\t"+str(p300),sep="\t",flush=True)
reader.close()
os.remove(tempFile)

