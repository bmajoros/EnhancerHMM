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
import math
from BedReader import BedReader
from Pipe import Pipe
import SumLogProbs

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
    cmd=MUMMIE+"/fastb-slice "+SCHEMA+" "+fastbDir+"/"+record.chr+"."+timepoint+".fastb "+str(record.interval.begin)+" "+str(record.interval.end)+" "+outfile
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

nextID=1
reader=BedReader(P300_BED)
while(True):
    record=reader.nextRecord()
    if(not record): break
    sliceFastb(record,"t00","slice0.fastb")
    sliceFastb(record,"t3","slice3.fastb")
    posterior0=getPosterior("slice0.fastb")
    posterior3=getPosterior("slice3.fastb")
    #if(posterior0<0.1 and posterior3>0.9):
    if(posterior3>0.9):
        id="elem"+str(nextID)
        nextID+=1
        print(record.chr+"\t"+str(record.interval.begin)+"\t"+str(record.interval.end)+"\t"+id+"\t"+str(posterior3-posterior0),sep="\t")
reader.close()


