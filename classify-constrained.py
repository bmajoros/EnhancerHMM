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
from Rex import Rex
rex=Rex()

if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <pos-dir> <neg-dir> <pos.hmm> <neg.hmm> <min-len>\n")
(posDir,negDir,posHMM,negHMM,minLen)=sys.argv[1:]

#timepoint="t00"
#minLen=200
minHump=1
minPeak=1
#HALF_WIDTH=250
fgPrior=math.log(0.5)
bgPrior=math.log(0.5)
BASE="/home/bmajoros/GGR/bacs"
fgHMM=BASE+"/hmm/trained-pos5.hmm"
bgHMM=BASE+"/hmm/trained-neg1.hmm"
P300_BED=BASE+"/p300-bacs-white.bed"
MUMMIE=os.environ["MUMMIE"]
SCHEMA=BASE+"/hmm/p300.schema"
fastbDir=BASE+"/fastb"

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

def process(dir,posHMM,negHMM,label):
  files=os.listdir()
  n=len(files)
  for i in range(n):
    file=files[i]
    file=file.rstrip()
    if(not rex.find(".fastb",file)): continue
    numer=getLL(file,posHMM)
    denom=getLL(file,negHMM)
    ratio=numer-denom
    ROC.write(str(ratio)+"\t"+str(label)+"\n")

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
    if(not applyConstraints(fgHMM,tempFile,minHump,minPeak,minLen)):
        continue
    posterior=getPosterior(tempFile)
    if(posterior>0.9):
        id="elem"+str(nextID)
        nextID+=1
        print(record.chr+"\t"+str(record.interval.begin)+"\t"+str(record.interval.end)+"\t"+id+"\t"+str(posterior),sep="\t",flush=True)
reader.close()
os.remove(tempFile)

open(ROC,">$ROC") || die "can't write to file: $ROC\n";
process($posDir,$posHMM,$negHMM,1);
process($negDir,$posHMM,$negHMM,0);
close(ROC);

#########################################################################
